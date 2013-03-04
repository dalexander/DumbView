
class Consensus(object):
    def __init__(self, refWindow, sequence, confidence, coverage):
        assert (len(sequence) ==
                len(confidence) ==
                len(coverage))
        self.refWindow  = refWindow
        self.sequence   = sequence
        self.confidence = confidence
        self.coverage   = coverage

    def maskSequence(referenceSequence, noEvidenceConsensusCall):
        assert len(referenceSequence) == len(self.sequence)
        pass

def variantsFromConsensus(refWindow, referenceSequence, consensusObject):
    pass
