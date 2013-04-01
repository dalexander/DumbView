from GenomicConsensus.quiver import quiver as q
from GenomicConsensus.quiver import utils  as qu
from GenomicConsensus.consensus import *
import ConsensusCore as cc
from .Window import *

K = 3
quiverConfig = q.QuiverConfig()
overlap = 5

def enlargedReferenceWindow(refWin, refContig, overlap):
    refId, refStart, refEnd = refWin
    contigLength = len(refContig)
    return Window(refId,
                  max(0, refStart - overlap),
                  min(refEnd + overlap + 1, contigLength))

def consensus(cmpH5, refWindow, rowNumbers, referenceTable):
    # identify the enlarged interval [-5, +5]
    refContig = referenceTable.byKey(refWindow.refId).sequence
    eWindow = enlargedReferenceWindow(refWindow, refContig, overlap)
    refSeqInEnlargedWindow = refContig[eWindow.start:eWindow.end]

    # find 3-spanned intervals in the enlarged interval
    # call css for each interval
    subConsensi = []
    tStart = cmpH5.tStart[rowNumbers]
    tEnd   = cmpH5.tEnd[rowNumbers]
    coveredIntervals = qu.kSpannedIntervals(eWindow, K, tStart, tEnd)
    holes = qu.holes(eWindow, coveredIntervals)

    print coveredIntervals
    print holes

    for interval in sorted(coveredIntervals + holes):
        subWin = qu.subWindow(eWindow, interval)
        intStart, intEnd = interval
        intRefSeq = refContig[intStart:intEnd]

        if interval in coveredIntervals:
            clippedAlns = [ aln.clippedTo(*interval)
                            for aln in cmpH5[rowNumbers] ]
            css_ = qu.quiverConsensusForAlignments(subWin,
                                                  intRefSeq,
                                                  clippedAlns,
                                                  quiverConfig)
        else:
            css_ = noCallAsConsensus(subWin, intRefSeq)

        subConsensi.append(css_)

    # join subconsensus objects
    css = join(subConsensi)

    # align css back to refWindow, and clip
    ga = cc.Align(refSeqInEnlargedWindow, css_.sequence)
    targetPositions = cc.TargetToQueryPositions(ga)
    cssStart = targetPositions[refWindow.start-eWindow.start]
    cssEnd   = targetPositions[refWindow.end-eWindow.start]

    cssSequence    = css_.sequence[cssStart:cssEnd]
    cssQv          = css_.confidence[cssStart:cssEnd]

    consensusObj = Consensus(refWindow,
                             cssSequence,
                             cssQv)

    return consensusObj
