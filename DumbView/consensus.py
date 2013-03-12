import numpy as np

class Consensus(object):
    def __init__(self, refWindow, sequence, confidence, coverage):
        assert (len(sequence) ==
                len(confidence) ==
                len(coverage))
        self.refWindow  = refWindow
        self.sequence   = sequence
        self.confidence = confidence
        self.coverage   = coverage


#
# Functions that produce a consensus sequence string for a window
# where there is inadequate evidence.  Each takes as arguments a
# reference window, and the sequence of the entire associated
# reference contig.
#

def noCallAsConsensus(referenceSequence):
    return "N" * len(referenceSequence)

def referenceAsConsensus(referenceSequence):
    return referenceSequence.upper()

def lowercaseReferenceAsConsensus(referenceSequence):
    return referenceSequence.lower()

noEvidenceConsensusFactoryByName = \
    { "nocall"             : noCallAsConsensus,
      "reference"          : referenceAsConsensus,
      "lowercasereference" : lowercaseReferenceAsConsensus}


#
# Variant finding.  Do this once in GenomicConsensus
#

def variantsFromConsensus(refWindow, referenceSequence, consensusObject):
    pass
