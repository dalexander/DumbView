#
# A tool for viewing the alignments in a given reference window.
# Will be migrated to pbh5tools when possible.
#

import argparse, os, shlex, sys
from pbcore.io import CmpH5Reader, GffReader, FastaTable
from pbcore.util.ToolRunner import PBToolRunner
from DumbView.format import *
from DumbView.Window import *
from GenomicConsensus.utils import readsInWindow

def loadReferences(fastaFilename, cmpH5):
    return FastaTable(fastaFilename)

def extractCmpH5AndReferenceFromGff(gffReader):
    #
    # New way
    #
    cmpH5 = None
    reference = None
    for h in gffReader.headers:
        if h.startswith("##source-alignment-file"):
            cmpH5 = h.split()[1]
        elif h.startswith("##source-reference-file"):
            reference = h.split()[1]
        if cmpH5 and reference:
            return cmpH5, reference

    #
    # Old way
    #
    # This code is a horrible hack and an affront to good taste and I
    # blame Python's stdlib for making this the way I had to do it.
    cmpH5 = None
    reference = None
    for header in gffReader.headers:
        if header.startswith("##source-commandline"):
            args = shlex.split(header)
            for flag in ["-r", "--referenceFilename"]:
                if flag in args:
                    reference = args[args.index(flag) + 1]
                    break
                else:
                    for arg in args:
                        if arg.startswith(flag):
                            reference = arg.split("=")[1]
                            break
            for arg in args:
                if arg.endswith(".cmp.h5"):
                    cmpH5 = arg
                    break
    return cmpH5, reference

def mainGff(options):
    reader = GffReader(options.inputGff)
    cmpH5Fname, referenceFname = extractCmpH5AndReferenceFromGff(reader)
    # Allow overriding
    cmpH5Fname = options.inputCmpH5 or cmpH5Fname
    referenceFname = options.referenceFilename or referenceFname

    assert cmpH5Fname
    assert referenceFname

    cmpH5 = CmpH5Reader(cmpH5Fname)
    referenceTable = loadReferences(referenceFname, cmpH5)

    for gffRecord in reader:
        referenceSeq = gffRecord.get("reference", "-")
        variantSeq   = gffRecord.get("variantSeq", "-")
        variantConfidence = gffRecord.confidence
        variantSummary = "(%s > %s)" % (referenceSeq, variantSeq)
        print gffRecord.type, gffRecord.seqid, gffRecord.start, gffRecord.end, \
            variantSummary, variantConfidence
        refId = cmpH5.referenceInfo(gffRecord.seqid).ID
        refLength = cmpH5.referenceInfo(gffRecord.seqid).Length
        refWindow = makeDisplayWindow(refLength, options.width,
                                       Window(refId,
                                              gffRecord.start-10,
                                              gffRecord.end+10))
        if "rows" in gffRecord.attributes:
            rowNumbers = map(int, gffRecord.rows.split(","))
        else:
            rowNumbers = readsInWindow(cmpH5, refWindow, options.depth,
                                       minMapQV=options.minMapQV, strategy=options.sorting)
        formatWindow(cmpH5, refWindow, rowNumbers, referenceTable,
                     aligned=(gffRecord.type != "insertion"),
                     consensus=options.consensus, useColor=options.color)
        print

def mainCmpH5(options):
    cmpH5 = CmpH5Reader(options.inputCmpH5)
    if options.referenceFilename:
        referenceTable = loadReferences(options.referenceFilename, cmpH5)
    else:
        referenceTable = None

    for refWindow in options.referenceWindows:
        refId = cmpH5.referenceInfo(refWindow.refId).ID
        refName = cmpH5.referenceInfo(refWindow.refId).FullName
        refLength = cmpH5.referenceInfo(refWindow.refId).Length
        refWindow = refWindow._replace(refId=refId)
        refWindow = makeDisplayWindow(refLength, options.width, refWindow)

        if options.rowNumbers != None:
            rowNumbers = options.rowNumbers
        else:
            rowNumbers = readsInWindow(cmpH5, refWindow, options.depth,
                                       minMapQV=options.minMapQV, strategy=options.sorting)

        print windowToGffString(Window(refName, refWindow.start, refWindow.end))

        if options.oneAtATime:
            formatIndividualAlignments(cmpH5, refWindow, rowNumbers)
        else:
            formatWindow(cmpH5, refWindow, rowNumbers,
                         referenceTable, options.aligned, options.color,
                         options.consensus)
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

        self.parser.set_defaults(consensus=True)
        arg("--consensus",   "-C", action="store_true", dest="consensus")
        arg("--noConsensus", "-N", action="store_false", dest="consensus")

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
            elif fname.endswith(".cmp.h5"):
                self.args.inputCmpH5 = fname
            else:
                print "Invalid input file"
                sys.exit(-1)


    def getVersion(self):
        return "0.2"

    def run(self):
        _main(self.args)
        return 0
