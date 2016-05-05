#
# A tool for viewing the alignments in a given reference window.
# Will be migrated to pbh5tools when possible.
#

import argparse, os, os.path, shlex, sys
from pbcore.io import AlignmentSet, ReferenceSet, GffReader, IndexedFastaReader
from pbcore.util.ToolRunner import PBToolRunner
from DumbView.format import *
from DumbView.window import *
from DumbView import __VERSION__
from GenomicConsensus.utils import readsInWindow

def loadReferences(fastaFilename, alnReader):
    # as of 3.0, quiver can be called with a "ReferenceSet" XML
    # instead of just a FASTA.  Let's just unwrap the underlying FASTA
    # file.  This code still works if a FASTA was provided.
    dset = ReferenceSet(fastaFilename)
    fastas = dset.toExternalFiles()
    assert len(fastas) == 1
    return IndexedFastaReader(fastas[0])


def dumpVariantCsv(fname, alnReader, alns, gffRecord, width=5):
    # Dump a CSV file for pulserecognizer with the following columns:
    # Movie,HoleNumber,rStart,rEnd
    refId    = gffRecord.seqid
    refStart = gffRecord.start - 1  # 1 to 0 based coordinates
    refEnd   = gffRecord.end

    expandedRange = (refStart - width, refEnd + width)

    with open(fname, "w") as f:
        f.write("Movie,HoleNumber,rStart,rEnd\n")
        for aln in alns:
            if aln.spansReferenceRange(*expandedRange):
                ca = aln.clippedTo(*expandedRange)
                f.write("%s,%d,%d,%d\n" %
                        (ca.movieInfo.Name,
                         ca.HoleNumber,
                         ca.rStart,
                         ca.rEnd))


def formatVariantCsvLink(fname):
    print "  Link for pulse recognizer: <a href=./" + fname + ">" + fname + "</a>\n"

def extractCmpH5AndReferenceFromGff(gffReader):
    #
    # New way
    #
    alnReader = None
    reference = None
    for h in gffReader.headers:
        if h.startswith("##source-alignment-file"):
            alnReader = h.split()[1]
        elif h.startswith("##source-reference-file"):
            reference = h.split()[1]
        if alnReader and reference:
            return alnReader, reference

def mainGff(options):
    reader = GffReader(options.inputGff)
    alnsFname, referenceFname = extractCmpH5AndReferenceFromGff(reader)
    # Allow overriding
    alnsFname = options.inputCmpH5 or alnsFname
    referenceFname = options.referenceFilename or referenceFname

    assert os.path.isfile(alnsFname)
    assert os.path.isfile(referenceFname)

    alnReader = AlignmentSet(alnsFname, referenceFastaFname=referenceFname)

    if options.fofn is not None:
        alnReader.attach(options.fofn)

    referenceTable = loadReferences(referenceFname, alnReader)

    for i, gffRecord in enumerate(reader):
        referenceSeq = gffRecord.get("reference", "-")
        variantSeq   = gffRecord.get("variantSeq", "-")
        variantConfidence = gffRecord.confidence
        variantSummary = "(%s > %s)" % (referenceSeq, variantSeq)
        print gffRecord.type, gffRecord.seqid, gffRecord.start, gffRecord.end, \
            variantSummary, variantConfidence
        refId = gffRecord.seqid
        refLength = alnReader.referenceInfo(gffRecord.seqid).Length
        refWindow = makeDisplayWindow(refLength, options.width,
                                       Window(refId,
                                              gffRecord.start-10,
                                              gffRecord.end+10))
        if "rows" in gffRecord.attributes:
            alns = alnReader[map(int, gffRecord.rows.split(","))]
        else:
            alns = readsInWindow(alnReader, refWindow, options.depth,
                                 minMapQV=options.minMapQV, strategy=options.sorting)
        formatWindow(alnReader, refWindow, alns, referenceTable,
                     aligned=(gffRecord.type != "insertion"),
                     consensus=options.consensus,
                     useColor=options.color,
                     realign=options.realign)

        if options.pulseRecognizer:
            # CSV output for pulse recognizer
            print
            csvFname = "variant-" + str(i) +  ".csv"
            dumpVariantCsv(csvFname, alnReader, alns, gffRecord)
            formatVariantCsvLink(csvFname)

        print

def mainCmpH5(options):
    alnReader = AlignmentSet(options.inputCmpH5,
                             referenceFastaFname=options.referenceFilename)
    if options.fofn is not None:
        alnReader.attach(options.fofn)

    if options.referenceFilename:
        referenceTable = loadReferences(options.referenceFilename, alnReader)
    else:
        referenceTable = None

    for refWindow in options.referenceWindows:
        refId = refWindow.refId
        refName = alnReader.referenceInfo(refWindow.refId).FullName
        refLength = alnReader.referenceInfo(refWindow.refId).Length
        refWindow = refWindow._replace(refId=refId)
        refWindow = makeDisplayWindow(refLength, options.width, refWindow)

        if options.rowNumbers != None:
            alns = alnReader[options.rowNumbers]
        else:
            alns = readsInWindow(alnReader, refWindow, options.depth,
                                       minMapQV=options.minMapQV, strategy=options.sorting)

        print windowToGffString(Window(refName, refWindow.start, refWindow.end))

        if options.oneAtATime:
            formatIndividualAlignments(alnReader, refWindow, alns)
        else:
            formatWindow(alnReader, refWindow, alns,
                         referenceTable, options.aligned, options.color,
                         options.realign, options.consensus)
        print

def _main(options):
    if any([fn.endswith(".gff") or fn.endswith(".gff.gz")
            for fn in options.inputFilenames]):
        mainGff(options)
    else:
        mainCmpH5(options)

class DumbViewApp(PBToolRunner):

    def __init__(self):
        desc = "Command-line PacBio genome browser"
        super(DumbViewApp, self).__init__(desc)

        arg = self.parser.add_argument
        arg("inputFilenames", nargs="+", type=str, help=".cmp.h5 or .gff filename, or both")
        arg("--referenceWindows", "-w", type=windowsFromGffStrings, default=[])
        arg("--referenceFilename", "-r", default=None)
        arg("--depth", "-D", "-X", type=int, default=20)
        arg("--width", "-W", type=int, default=40)
        arg("--minMapQV", "-m", type=int, default=10)
        arg("--rowNumbers", "-n", type=int, nargs="+", default=None)
        arg("--columns", type=str, nargs="+", default=None)
        arg("--unaligned", "-u", dest="aligned", action="store_false")
        arg("--aligned",   "-a", dest="aligned", action="store_true", default=True)
        arg("--oneAtATime", "-1", action="store_true", default=False)
        arg("--sorting", "-s", choices=["fileorder", "longest", "spanning"], default="longest")
        arg("--fofn", default=None)
        arg("--pulseRecognizer", action="store_true", default=False,
            help="In variants.gff analysis mode, emit a CSV file for inspection in PulseRecognizer")
        self.parser.set_defaults(consensus=True)
        arg("--consensus",   "-C", action="store_true", dest="consensus")
        arg("--noConsensus", "-N", action="store_false", dest="consensus")
        self.parser.set_defaults(realign=True)
        arg("--realign",   action="store_true", dest="realign", help="Enable simple gap-pushing realignment")
        arg("--noRealign", action="store_false", dest="realign", help="Disable gap-pushing realignment; alignments directly from file")

        class ColorAction(argparse.Action):
            def __call__(self, parser, namespace, values, option_string=None):
                if (values == None) or (values == "on"):
                    color = True
                elif (values == "off"):
                    color = False
                else:
                    assert values == "auto"
                    color = os.isatty(1)
                setattr(namespace, self.dest, color)

        arg("--color", "-c", nargs="?", choices=["on", "off", "auto"],
            action=ColorAction, default=os.isatty(1))


    def validateArgs(self):
        self.args.inputGff  = None
        self.args.inputCmpH5 = None
        for fname in self.args.inputFilenames:
            if fname.endswith(".gff") or fname.endswith(".gff.gz"):
                self.args.inputGff = fname
            elif fname.endswith(".cmp.h5") or fname.endswith(".bam") or fname.endswith(".fofn"):
                self.args.inputCmpH5 = fname
            else:
                print "Invalid input file"
                sys.exit(-1)


    def getVersion(self):
        return __VERSION__

    def run(self):
        try:
            import ipdb
            with ipdb.launch_ipdb_on_exception():
                _main(self.args)
            return 0
        except ImportError:
            _main(self.args)
