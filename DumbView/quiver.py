
import numpy as np
import ConsensusCore as cc
from GenomicConsensus.quiver import model, utils as qutils

from DumbView.consensus import Consensus, noCallAsConsensus
from DumbView.utils import readsInWindow, kSpannedIntervals, holes


class QuiverConfig(object):
    def __init__(self,
                 minMapQV=10,
                 minPoaCoverage=3,
                 maxPoaCoverage=11,
                 mutationSeparation=10,
                 maxIterations=20,
                 refineDinucleotideRepeats=True,
                 noEvidenceConsensus=noCallAsConsensus,
                 parameters=None):

        self.minMapQV                   = minMapQV
        self.minPoaCoverage             = minPoaCoverage
        self.maxPoaCoverage             = maxPoaCoverage
        self.mutationSeparation         = mutationSeparation
        self.maxIterations              = maxIterations
        self.refineDinucleotideRepeats  = refineDinucleotideRepeats
        self.noEvidenceConsensus        = noEvidenceConsensus
        self.parameters                 = parameters

        # Convenience
        self.model                      = self.parameters.model
        self.ccQuiverConfig             = self.parameters.quiverConfig


    @staticmethod
    def fromOptions(options, cmpH5):
        pass

def quiverConsensusForWindow(cmpH5, refWindow, referenceContig,
                             depthLimit, quiverConfig):
    """
    High-level routine for calling the consensus for a
    window of the genome given a cmp.h5.

    Identifies the coverage contours of the window in order to
    identify subintervals where a good consensus can be called.
    Creates the desired "no evidence consensus" where there is
    inadequate coverage.
    """
    # 1) identify the intervals with adequate coverage for quiver
    #    consensus; restrict to intervals of length > 10
    allRows = readsInWindow(cmpH5, refWindow, minMapQV=quiverConfig.minMapQV)
    starts = cmpH5.tStart[allRows]
    ends   = cmpH5.tEnd[allRows]
    intervals = [ (s, e)
                  for (s, e) in kSpannedIntervals(refWindow,
                                                  quiverConfig.minPoaCoverage,
                                                  starts,
                                                  ends)
                  if (e - s) > 10 ]

    coverageGaps = holes(refWindow, intervals)
    allIntervals = sorted(intervals + coverageGaps)
    assert holes(refWindow, allIntervals) == []

    # 2) pull out the reads we will use for each interval
    # 3) call quiverConsensusForRows on the interval
    subConsensi = []
    for interval in allIntervals:
        intStart, intEnd = interval
        intRefSeq = referenceContig[intStart:intEnd]
        subWindow = refWindow.subWindow(*interval)

        if interval in coverageGaps:
            cssSeq = quiverConfig.noEvidenceConsensus(intRefSeq)
            css = Consensus(subWindow,
                            cssSeq,
                            [0]*len(cssSeq),
                            [0]*len(cssSeq))
        else:

            windowRefSeq = referenceContig[intStart:intEnd]
            rows = readsInWindow(cmpH5, subWindow,
                                 depthLimit=depthLimit,
                                 minMapQV=quiverConfig.minMapQV,
                                 strategy="longest")

            # TODO: Some further filtering: remove "stumpy reads"
            alns = cmpH5[rows]
            clippedAlns = [ aln.clippedTo(*interval) for aln in alns ]
            css = quiverConsensusForAlignments(subWindow,
                                               intRefSeq,
                                               clippedAlns,
                                               quiverConfig)
        subConsensi.append(css)

    # 4) glue the subwindow consensus objects together to form the
    #    full window consensus
    cssSeq_  = "".join(sc.sequence for sc in subConsensi)
    cssConf_ = np.concatenate([sc.confidence for sc in subConsensi])
    cssCov_  = np.concatenate([sc.coverage for sc in subConsensi])
    return Consensus(refWindow,
                     cssSeq_,
                     cssConf_,
                     cssCov_)


def quiverConsensusForAlignments(refWindow, refSequence, alns, quiverConfig):
    """
    Testable!

    refSequence is reference within window.
    Q: is it needed for anything beyond reporting?
       if not, we should remove it from prototype

    Clipping has already been done!
    """
    # (The buck stops here.  Call consensus on this interval, no recursing)
    refStart = refWindow.start
    refEnd   = refWindow.end

    # Compute the POA consensus, which is our initial guess, and
    # should typically be > 99.5% accurate
    fwdSequences = [ a.read(orientation="genomic", aligned=False)
                     for a in alns
                     if a.spansReferenceRange(refStart, refEnd) ]
    assert len(fwdSequences) >= quiverConfig.minPoaCoverage
    p = cc.PoaConsensus.FindConsensus(fwdSequences[:quiverConfig.maxPoaCoverage])
    ga = cc.Align(refSequence, p.Sequence())
    numPoaVariants = ga.Errors()
    poaCss = p.Sequence()

    # Extract reads into ConsensusCore-compatible objects, and map them into the
    # coordinates relative to the POA consensus
    mappedReads = [ quiverConfig.model.extractMappedRead(aln, refStart)
                    for aln in alns ]
    queryPositions = cc.TargetToQueryPositions(ga)
    mappedReads = [ qutils.lifted(queryPositions, mr) for mr in mappedReads ]

    # Load the mapped reads into the mutation scorer, and iterate
    # until convergence.
    mms = cc.SparseSseQvMultiReadMutationScorer(quiverConfig.ccQuiverConfig, poaCss)
    for mr in mappedReads:
        mms.AddRead(mr)

    # Iterate until covergence
    # TODO: pass quiverConfig down here.
    _, quiverConverged = qutils.refineConsensus(mms)
    if quiverConfig.refineDinucleotideRepeats:
        qutils.refineDinucleotideRepeats(mms)
    quiverCss = mms.Template()

    # TODO: calculate confidence here
    confidence = [0] * len(quiverCss)
    coverage   = [0] * len(quiverCss)

    return Consensus(refWindow,
                     quiverCss,
                     confidence,
                     coverage)
