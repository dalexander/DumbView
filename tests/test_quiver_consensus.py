
from pbcore.io import FastaTable, CmpH5Reader

from DumbView.quiver import *
from DumbView.utils import Window
from GenomicConsensus.quiver.model import *
from AlignmentHitStubs import ForwardAndReverseReads



def test_quiver_1():
    p = loadParameterSets(findParametersFile())["unknown.AllQVsModel"]
    myQuiverConfig = QuiverConfig(parameters=p)
    css = quiverConsensusForAlignments(Window(0, 0, 20),
                                       ForwardAndReverseReads.reference,
                                       ForwardAndReverseReads.hits,
                                       myQuiverConfig)
    assert css.sequence == ForwardAndReverseReads.reference[0:20]


def test_quiver_2():
    p = loadParameterSets(findParametersFile())["unknown.AllQVsModel"]
    myQuiverConfig = QuiverConfig(parameters=p)

    c = CmpH5Reader("~/Data/038537.cmp.h5")
    fastaTable = FastaTable("~/Data/lambdaNEB.fa")
    seq = fastaTable[0].sequence[0:60]

    #alns = c[0:12] + c[13:30]
    css = quiverConsensusForWindow(c, Window("ref000001", 0, 60), seq,
                                   20, myQuiverConfig)
