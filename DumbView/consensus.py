import GenomicConsensus.quiver as q
from GenomicConsensus.consensus import *
import ConsensusCore as cc
from .Window import *
from .utils import readsInWindow

K = 3
quiverConfig = q.model.QuiverConfig()
overlap = 5

def enlargedReferenceWindow(refWin, contigLength, overlap):
    refId, refStart, refEnd = refWin
    return Window(refId,
                  max(0, refStart - overlap),
                  min(refEnd + overlap + 1, contigLength))

def consensus(cmpH5, refWindow, referenceTable, rowNumbers=None):
    # identify the enlarged interval [-5, +5]
    refName = cmpH5.referenceInfo(refWindow.refId).FullName
    refLength = referenceTable.length(refName)
    eWindow = enlargedReferenceWindow(refWindow, refLength, overlap)
    refSeqInWindow = referenceTable.sequence(refName, refWindow.start, refWindow.end)
    refSeqInEnlargedWindow = referenceTable.sequence(refName, eWindow.start, eWindow.end)

    # find 3-spanned intervals in the enlarged interval
    # call css for each interval
    subConsensi = []
    tStart = cmpH5.tStart[rowNumbers]
    tEnd   = cmpH5.tEnd[rowNumbers]
    coveredIntervals = q.utils.kSpannedIntervals(eWindow, K, tStart, tEnd)
    holes = q.utils.holes(eWindow, coveredIntervals)

    for interval in sorted(coveredIntervals + holes):
        subWin = subWindow(eWindow, interval)
        #print subWin
        intStart, intEnd = interval
        intRefSeq = refSeqInEnlargedWindow[intStart-eWindow.start:
                                           intEnd-eWindow.start]
        css_ = noCallAsConsensus(subWin, intRefSeq)
        if interval in coveredIntervals:
            rows = readsInWindow(cmpH5, subWin,
                                 depthLimit=100,
                                 minMapQV=quiverConfig.minMapQV,
                                 strategy="longest")
            clippedAlns = [ aln.clippedTo(*interval) for aln in cmpH5[rows]]
            goodAlns = q.utils.filterAlnsForQuiver(subWin, clippedAlns, quiverConfig)
            if len(goodAlns) >= K:
                css_ = q.utils.quiverConsensusForAlignments(subWin,
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
    ga = cc.AlignWithAffineGapPenalty(ref, query)
    return (ga.Target(),
            ga.Transcript(),
            ga.Query())
