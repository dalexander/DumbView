
from DumbView.consensus import Consensus

import numpy as np
from itertools import izip
from collections import Counter


def pluralityConsensusForRows(cmpH5, refWindow,
                              rowNumbers, realignHomopolymers=False):
    alns = cmpH5[rowNumbers]
    return pluralityConsensusFromAlignments(refWindow, alns, realignHomopolymers)

def pluralityConsensusFromAlignments(refWindow, alns, realignHomopolymers=False):
    """
    Compute a Consensus object for this window, using the given
    `alns`, by applying a straightforward column-oriented consensus
    calling algorithm.

    If the consensus cannot be called for a base, "N" will be placed
    in the consensus sequence for that position.

    If `realignHomopolymers` is True, alignment gaps will be shuffled
    in homopolymer regions in an attempt to maximize variant detection
    sensitivity.
    """
    refStart = refWindow.start
    refEnd = refWindow.end
    windowSize = refEnd - refStart

    baseCallsMatrix = tabulateBaseCalls(refWindow, alns, realignHomopolymers)

    #
    # Summarize the baseCallsMatrix table, creating a Consensus object
    #
    consensusSequence  = np.empty(shape=windowSize, dtype="S1")
    consensusFrequency = np.empty(shape=windowSize, dtype=int)
    effectiveCoverage  = np.empty(shape=windowSize, dtype=int)
    for j in xrange(0, windowSize):
        counter = Counter(baseCallsMatrix[:, j])
        if "" in counter: counter.pop("")
        effectiveCoverage[j] = sum(counter.itervalues())
        if effectiveCoverage[j] == 0:
            consensusFrequency[j], consensusSequence[j] = 0, "N"
        else:
            consensusFrequency[j], consensusSequence[j] = max(zip(counter.values(),
                                                                  counter.keys()))
    consensusConfidence = computePluralityConfidence(consensusFrequency,
                                                     effectiveCoverage)
    return PluralityConsensus(refWindow,
                              consensusSequence.tostring(),
                              consensusConfidence,
                              effectiveCoverage,
                              consensusFrequency)



class PluralityConsensus(Consensus):
    def __init__(self, refWindow, sequence, confidence, coverage, frequency):
        assert len(frequency) == len(sequence)
        super(PluralityConsensus, self).__init__(refWindow,
                                                 sequence,
                                                 confidence,
                                                 coverage)
        self.frequency = frequency


def tabulateBaseCalls(refWindow, alns, realignHomopolymers=False):
    #
    # Go through the reads and build up the structured baseCallsMatrix table
    #
    refStart = refWindow.start
    refEnd = refWindow.end
    windowSize = refEnd - refStart

    baseCallsMatrix = np.zeros(shape=(len(alns), windowSize), dtype="S8")

    for i, aln in enumerate(alns):
        aln = aln.clippedTo(refStart, refEnd)
        alnRef    = aln.reference(orientation="genomic")
        alnRead   = aln.read(orientation="genomic")
        alnRefPos = aln.referencePositions(orientation="genomic")
        if realignHomopolymers:
            alnRef, alnRead, alnRefPos = \
                normalizeHomopolymerGaps(alnRef, alnRead, alnRefPos)

        readBases = []
        for (refPos, refBase, readBase) in izip(alnRefPos,
                                                alnRef,
                                                alnRead):
            if readBase != "-":
                readBases.append(readBase)

            if refBase != "-":
                if not readBases:
                    basesForRefPos = "-"
                else:
                    basesForRefPos = "".join(readBases)
                baseCallsMatrix[i, (refPos - refStart)] = basesForRefPos
                readBases = []

    return baseCallsMatrix

def computePluralityConfidence(consensusFrequency, effectiveCoverage):
    """
    Come up with a new simpler scheme.  Should only need frequency and
    coverage.
    """
    confidence = np.empty_like(consensusFrequency)
    confidence.fill(25)
    return confidence

def normalizeHomopolymerGaps(alnRef, alnRead, alnRefPos=None):
    pass
