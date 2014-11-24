import GenomicConsensus.quiver as q
import ConsensusCore as cc
import GenomicConsensus.windows as w

from GenomicConsensus.consensus import *
from GenomicConsensus.utils import readsInWindow

from .window import Window, subWindow


K = 3
quiverConfig = q.model.loadQuiverConfig("unknown.NoQVsModel")
overlap = 5

def enlargedReferenceWindow(refWin, contigLength, overlap):
    refId, refStart, refEnd = refWin
    return Window(refId,
                  max(0, refStart - overlap),
                  min(refEnd + overlap + 1, contigLength))

def consensus(cmpH5, refWindow, referenceTable, alns):
    # identify the enlarged interval [-5, +5]
    refName = cmpH5.referenceInfo(refWindow.refId).FullName
    refLength = len(referenceTable[refName].sequence)
    eWindow = enlargedReferenceWindow(refWindow, refLength, overlap)
    refSeqInEnlargedWindow = referenceTable[refName].sequence[eWindow.start:eWindow.end]

    # find 3-spanned intervals in the enlarged interval
    # call css for each interval
    subConsensi = []
    tStart = [ a.tStart for a in alns ]
    tEnd = [ a.tEnd for a in alns ]
    coveredIntervals = w.kSpannedIntervals(eWindow, K, tStart, tEnd)
    holes = w.holes(eWindow, coveredIntervals)

    for interval in sorted(coveredIntervals + holes):
        subWin = subWindow(eWindow, interval)
        #print subWin
        intStart, intEnd = interval
        intRefSeq = refSeqInEnlargedWindow[intStart-eWindow.start:
                                           intEnd-eWindow.start]
        css_ = Consensus.nAsConsensus(subWin, intRefSeq)
        if interval in coveredIntervals:
            alns = readsInWindow(cmpH5, subWin,
                                 depthLimit=100,
                                 minMapQV=quiverConfig.minMapQV,
                                 strategy="longest")
            clippedAlns = [ aln.clippedTo(*interval) for aln in alns ]
            goodAlns = q.utils.filterAlns(subWin, clippedAlns, quiverConfig)
            if len(goodAlns) >= K:
                css_ = q.utils.consensusForAlignments(subWin,
                                                      intRefSeq,
                                                      goodAlns,
                                                      quiverConfig)

        subConsensi.append(css_)

    # join subconsensus objects
    css = join(subConsensi)

    # align css back to refWindow, and clip
    ga = cc.Align(refSeqInEnlargedWindow, css.sequence)
    targetPositions = cc.TargetToQueryPositions(ga)
    cssStart = targetPositions[refWindow.start-eWindow.start]
    cssEnd   = targetPositions[refWindow.end-eWindow.start]

    cssSequence    = css.sequence[cssStart:cssEnd]
    cssQv          = css.confidence[cssStart:cssEnd]

    consensusObj = Consensus(refWindow,
                             cssSequence,
                             cssQv)
    return consensusObj


def align(ref, query):
    ga = cc.AlignAffine(ref, query)
    return (ga.Target(),
            ga.Transcript(),
            ga.Query())
