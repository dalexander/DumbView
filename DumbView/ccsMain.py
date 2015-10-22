from pbcore.io import loadAlignmentFile, FastaTable
from pbcore.util.ToolRunner import PBToolRunner
import sys, argparse, os

from DumbView.main import loadReferences
from DumbView.window import Window
from DumbView.format import *
from DumbView import __VERSION__

def windowChunks(refWindow, WIDTH=60):
    refId = refWindow.refId
    for s in xrange(refWindow.start, refWindow.end, WIDTH):
        e = min(refWindow.end, s + WIDTH)
        yield Window(refId, s, e)

def main(options):
    alnSubreads = loadAlignmentFile(options.subreads)
    alnCCS = loadAlignmentFile(options.ccs)
    if options.referenceFilename:
        referenceTable = loadReferences(options.referenceFilename, alnSubreads)
    else:
        referenceTable = None

    subreads = alnSubreads.readsForZmw(options.zmw)
    ccsReads = alnCCS.readsForZmw(options.zmw)
    allReads = subreads + ccsReads

    assert len(ccsReads) == 1
    assert len(subreads) > 0

    ccsRead = ccsReads[0]

    refIds = set([a.referenceId for a in allReads])
    assert len(refIds) == 1
    refId = ccsRead.referenceId
    refStart = min(a.referenceStart for a in allReads)
    refEnd = max(a.referenceEnd for a in allReads)

    refWindow = Window(refId, refStart, refEnd)
    ccsRowNumber = ccsRead.rowNumber

    if options.oneAtATime:
        formatIndividualAlignments(alnSubreads, refWindow, subreads)

    else:
        for subWindow in windowChunks(refWindow):
            formatWindow(subreadCmp, subWindow, subreads,
                         referenceTable, options.aligned, options.color, options.realign)
            print


class DumbViewCCSApp(PBToolRunner):

    def __init__(self):
        desc = "Command-line PacBio genome browser"
        super(DumbViewCCSApp, self).__init__(desc)

        arg = self.parser.add_argument
        arg("--subreads", type=str, required=True)
        arg("--ccs", type=str, required=True)
        arg("--zmw", type=str, required=True)

        arg("--referenceFilename", "-r", default=None)
        arg("--width", "-W", type=int, default=40)
        arg("--minMapQV", "-m", type=int, default=10)
        arg("--rowNumbers", "-n", type=int, nargs="+", default=None)
        arg("--columns", type=str, nargs="+", default=None)
        arg("--unaligned", "-u", dest="aligned", action="store_false")
        arg("--aligned",   "-a", dest="aligned", action="store_true", default=True)
        arg("--oneAtATime", "-1", action="store_true", default=False)
        arg("--sorting", "-s", choices=["readorder", "strand", "longest"], default="readorder")

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

    def getVersion(self):
        return __VERSION__


    def run(self):
        try:
            import ipdb
            with ipdb.launch_ipdb_on_exception():
                main(self.args)
            return 0
        except ImportError:
            main(self.args)


if __name__ == "__main__":
    sys.exit(DumbViewCCSApp().start())
