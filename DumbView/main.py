#
# A tool for viewing the alignments in a given reference window.
# Will be migrated to pbh5tools when possible.
#

import argparse, cProfile, numpy as np, os, shlex, pstats, sys
from pbcore.io import CmpH5Reader, FastaReader, GffReader
from DumbView.format import *
from DumbView.utils import *
from DumbView.Window import *
from DumbView.FastaTable import *

def loadReferences(fastaFilename, cmpH5):
    return FastaTable(fastaFilename)

def parseOptions():
    parser = argparse.ArgumentParser(description="View alignments")
    parser.add_argument("inputFilename", type=str, help=".cmp.h5 or .gff filename")
    parser.add_argument("--referenceWindow", "-w", type=windowFromGffString, default=None)
    parser.add_argument("--referenceFilename", "-r", default=None)
    parser.add_argument("--depth", "-X", type=int, default=20)
    parser.add_argument("--minMapQV", "-m", type=int, default=10)
    parser.add_argument("--rowNumbers", type=int, nargs="+", default=None)
    parser.add_argument("--columns", type=str, nargs="+", default=None)
    parser.add_argument("--unaligned", "-u", dest="aligned", action="store_false")
    parser.add_argument("--aligned",   "-a", dest="aligned", action="store_true", default=True)
    parser.add_argument("--sorting", "-s", choices=["fileorder", "longest", "spanning"],
                        default="longest")
    parser.add_argument("--profile", action="store_true", dest="doProfiling")

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

    parser.add_argument("--color", "-c", nargs="?", choices=["on", "off", "auto"],
                        action=ColorAction, default=os.isatty(1))

    options = parser.parse_args()
    return options


def extractCmpH5AndReferenceFromGffHeader(header):
    # This code is a horrible hack and an affront to good taste and I
    # blame Python's stdlib for making this the way I had to do it.
    assert header.startswith("##source-commandline")
    cmpH5 = None
    reference = None
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
    reader = GffReader(options.inputFilename)
    for header in reader.headers:
        if header.startswith("##source-commandline"):
            cmpH5Fname, referenceFname = extractCmpH5AndReferenceFromGffHeader(header)

    print cmpH5Fname

    cmpH5 = CmpH5Reader(cmpH5Fname)
    referenceTable = loadReferences(referenceFname, cmpH5)

    for gffRecord in reader:
        referenceSeq = gffRecord.attributes.get("reference", "-")
        variantSeq   = gffRecord.attributes.get("variantSeq", "-")
        variantSummary = "(%s > %s)" % (referenceSeq, variantSeq)
        print gffRecord.type, gffRecord.seqid, gffRecord.start, gffRecord.end, variantSummary
        refId = cmpH5.referenceInfo(gffRecord.seqid).ID
        refWindow = Window(refId,
                           gffRecord.start - 10,
                           gffRecord.end   + 10)
        rowNumbers = readsInWindow(cmpH5, refWindow, options.depth,
                                   minMapQV=options.minMapQV, strategy=options.sorting)
        formatWindow(cmpH5, refWindow, rowNumbers, referenceTable,
                     aligned=(gffRecord.type != "insertion"))
        print

def mainCmpH5(options):
    cmpH5 = CmpH5Reader(options.inputFilename)
    # Canonicalize window here
    refId = cmpH5.referenceInfo(options.referenceWindow.refId).ID
    refWindow = options.referenceWindow

    if options.rowNumbers != None:
        rowNumbers = options.rowNumbers
    else:
        rowNumbers = readsInWindow(cmpH5, refWindow, options.depth,
                                   minMapQV=options.minMapQV, strategy=options.sorting)

    if options.referenceFilename:
        referenceTable = loadReferences(options.referenceFilename, cmpH5)
    else:
        referenceTable = None

    formatWindow(cmpH5, refWindow, rowNumbers,
                 referenceTable, options.aligned, options.color)

def _main(options):
    options = parseOptions()
    if options.inputFilename.endswith(".cmp.h5"):
        mainCmpH5(options)
    elif options.inputFilename.endswith(".gff"):
        mainGff(options)
    else:
        print "Invalid input filename!"
        sys.exit(1)

def main():
    options = parseOptions()
    if options.doProfiling:
        cProfile.runctx("_main(options)", globals(), locals(), "profile.out")
        pstats.Stats("profile.out").sort_stats("cumulative").print_stats(20)
    else:
        _main(options)
