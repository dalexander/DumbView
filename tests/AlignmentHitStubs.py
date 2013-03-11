"""
This module provides an `AlignmentHitStub` class as well as a group of
curated stub objects, allowing for decoupled testing.
"""

import numpy as np
from GenomicConsensus.utils import complement, reverseComplement
from GenomicConsensus.variants import *

def ungappedPulseArray(a):
    dtype = a.dtype
    if dtype == np.float32:
        return a[~np.isnan(a)]
    elif dtype == np.uint8:
        return a[a != np.uint8(-1)]
    elif dtype == np.uint16:
        return a[a != np.uint16(-1)]
    elif dtype == np.int8:
        return a[a != ord("-")]
    else:
        raise Exception, "Invalid pulse array type"

def _makePulseFeatureAccessor(featureName):
    def f(self, aligned=True, orientation="native"):
        return self.pulseFeature(featureName, aligned=aligned, orientation=orientation)
    return f

class AlignmentHitStub(object):

    def __init__(self, referenceStart, reverseStrand, nativeReference, read, **kwargs):
        """
        Initialize stub with data, as if it just came from the
        instrument (i.e., reference must be provided in
        reverse-complemented form for reverse-strand reads).
        """
        self.reverseStrand = reverseStrand
        self.forwardStrand = not reverseStrand
        self.RCRefStrand = 1 if reverseStrand else 0
        self.referenceStart = referenceStart
        self.referenceEnd = referenceStart + sum(b != '-' for b in nativeReference)
        self._reference = np.fromstring(nativeReference, dtype="S1")
        self._read      = np.fromstring(read, dtype="S1")

        self._pulseFeatures = {}
        for featureName, feature in kwargs.iteritems():
            self._pulseFeatures[featureName] = feature

    def read(self, aligned=True, orientation="native"):
        val = self._read
        if not aligned:
            val = val[val != "-"]
        if orientation == "genomic" and self.reverseStrand:
            val = reverseComplement(val)
        return val.tostring()

    def reference(self, aligned=True, orientation="native"):
        val = self._reference
        if not aligned:
            val = val[val != "-"]
        if orientation == "genomic" and self.reverseStrand:
            val = reverseComplement(val)
        return val.tostring()

    def referencePositions(self, orientation="native"):
        genomicReference = np.fromstring(self.reference(orientation="genomic"), dtype="S1")
        genomicPositions =        \
            self.referenceStart + \
            np.append(0, np.cumsum(genomicReference != "-")[:-1])
        if orientation == "native" and self.reverseStrand:
            return genomicPositions[::-1]
        else:
            return genomicPositions

    def spansReferencePosition(self, refPos):
        return self.referenceStart <= refPos < self.referenceEnd

    def spansReferenceRange(self, start, end):
        assert start <= end
        return (self.referenceStart <= start <= end <= self.referenceEnd)

    @property
    def referenceSpan(self):
        return self.referenceEnd - self.referenceStart

    def pulseFeature(self, featureName, aligned=True, orientation="native"):
        data = self._pulseFeatures[featureName]
        if orientation == "genomic" and self.reverseStrand:
            data = data[::-1]
        if not aligned:
            data = ungappedPulseArray(data)
        return data

    IPD                = _makePulseFeatureAccessor("IPD")
    PulseWidth         = _makePulseFeatureAccessor("PulseWidth")
    QualityValue       = _makePulseFeatureAccessor("QualityValue")
    InsertionQV        = _makePulseFeatureAccessor("InsertionQV")
    DeletionQV         = _makePulseFeatureAccessor("DeletionQV")
    DeletionTag        = _makePulseFeatureAccessor("DeletionTag")
    MergeQV            = _makePulseFeatureAccessor("MergeQV")
    SubstitutionQV     = _makePulseFeatureAccessor("SubstitutionQV")

    def __repr__(self):
        return "Stub 0x%x" % id(self)

    def __getattr__(self, key):
        pass


FORWARD, REVERSE = False, True

def _(s):
    """
    Decode an ASCII-art representation of
    a pulse feature.

    Spec: string -> np.array(dtype=float32)

        _(" ") -> [20],         error improbable
        _("*") -> [3],          error probable
        _("-") -> [NaN],        for an alignment gap
        _("A") -> [ord("A")],   (works for A,T,G,C,N)
    """
    _decoder = { " " : 20,
                 "*" : 0,
                 "-" : np.NaN }
    for c in "ATGCN":
        _decoder[c] = ord(c)
    return np.array([_decoder[c] for c in s], dtype=np.float32)




class ForwardAndReverseReads(object):
    """
    3 fwd, 3 rev; blunt ends.

    FWD READS:        ccaaaaccccc ttttggggcc  (hit1; has insertion)
                      ccaaaa cccc ttttg-ggcc  (hit2; has deletion)
                      ccagaa cccc ttttggggcc  (hit3; has substitution)
    REFERENCE:  (fwd) CCAAAA CCCC TTTTGGGGCC
                (pos) 012345 6789 0123456789
                (rev) GGTTTT GGGG AAAACCCCGG
    REV READS:        ggtttt gggggaaaaccccgg  (hit4 = hit1'; has ins)
                      ggtttt gggg aaaac-ccgg  (hit5 = hit2'; has del)
                      ggtctt gggg aaaaccccgg  (hit6 = hit3'; has sub)

    Notes:

    - hit4-hit6 sequences are the reverse complements of hit1-hit3;
      with gaps and insqv/delqv bumps placed where blasr+primary would
      put them.

    - No plurality variants.

    """
    reference                  = "CCAAAACCCCTTTTGGGGCC"
    expectedPluralityConsensus = "CCAAAACCCCTTTTGGGGCC"
    referenceWindow = (1, 0, 20)

    hit1 = AlignmentHitStub(0, FORWARD,
                            "CCAAAA-CCCCTTTTGGGGCC",
                            "CCAAAACCCCCTTTTGGGGCC",
          InsertionQV=    _("      *              "),
          SubstitutionQV= _("                     "),
          DeletionQV=     _("                     "),
          DeletionTag=    _("NNNNNNNNNNNNNNNNNNNNN"),
          MergeQV=        _("                     "))

    hit2 = AlignmentHitStub(0, FORWARD,
                            "CCAAAACCCCTTTTGGGGCC",
                            "CCAAAACCCCTTTTG-GGCC",
          InsertionQV=    _("               -    "),
          SubstitutionQV= _("               -    "),
          DeletionQV=     _("               -*   "),
          DeletionTag=    _("NNNNNNNNNNNNNNN-GNNN"),
          MergeQV=        _("               -    "))

    hit3 = AlignmentHitStub(0, FORWARD,
                            "CCAAAACCCCTTTTGGGGCC",
                            "CCAGAACCCCTTTTGGGGCC",
          InsertionQV=    _("                    "),
          SubstitutionQV= _("   *                "),
          DeletionQV=     _("                    "),
          DeletionTag=    _("NNNNNNNNNNNNNNNNNNNN"),
          MergeQV=        _("                    "))

    hit4 = AlignmentHitStub(0, REVERSE,
                            "GGTTTTGGGG-AAAACCCCGG"  [::-1],
                            "GGTTTTGGGGGAAAACCCCGG"  [::-1],
          InsertionQV=    _("          *          ") [::-1],
          SubstitutionQV= _("                     ") [::-1],
          DeletionQV=     _("                     ") [::-1],
          DeletionTag=    _("NNNNNNNNNNNNNNNNNNNNN") [::-1],
          MergeQV=        _("                     ") [::-1])

    hit5 = AlignmentHitStub(0, REVERSE,
                            "GGTTTTGGGGAAAACCCCGG"  [::-1],
                            "GGTTTTGGGGAAAAC-CCGG"  [::-1],
          InsertionQV=    _("               -    ") [::-1],
          SubstitutionQV= _("               -    ") [::-1],
          DeletionQV=     _("              *-    ") [::-1],
          DeletionTag=    _("NNNNNNNNNNNNNNC-NNNN") [::-1],
          MergeQV=        _("               -    ") [::-1])

    hit6 = AlignmentHitStub(0, REVERSE,
                            "GGTTTTGGGGAAAACCCCGG"  [::-1],
                            "GGTCTTGGGGAAAACCCCGG"  [::-1],
          InsertionQV=    _("                    ") [::-1],
          SubstitutionQV= _("   *                ") [::-1],
          DeletionQV=     _("                    ") [::-1],
          DeletionTag=    _("NNNNNNNNNNNNNNNNNNNN") [::-1],
          MergeQV=        _("                    ") [::-1])

    hits = [hit1, hit2, hit3, hit4, hit5, hit6]

    insertionVariant = Insertion(1, 6, 6, "-", "A")
    deletionVariant  = Deletion(1, 15, 16, "G", "-")
    substitutionVariant = Substitution(1, 3, 4, "A", "G")

    allVariants = set([insertionVariant,
                       deletionVariant,
                       substitutionVariant])
    pluralityVariants = set([])



class AllForwardStrandReads(object):
    """
    3 forward strand reads with blunt ends.
    """
    reference                  = ForwardAndReverseReads.reference
    expectedPluralityConsensus = ForwardAndReverseReads.expectedPluralityConsensus
    referenceWindow = ForwardAndReverseReads.referenceWindow
    hits = [ForwardAndReverseReads.hit1,
            ForwardAndReverseReads.hit2,
            ForwardAndReverseReads.hit3]

class AllReverseStrandReads(object):
    """
    3 reverse strand reads with blunt ends.
    """
    reference                  = ForwardAndReverseReads.reference
    expectedPluralityConsensus = ForwardAndReverseReads.expectedPluralityConsensus
    referenceWindow = ForwardAndReverseReads.referenceWindow
    hits = [ForwardAndReverseReads.hit4,
            ForwardAndReverseReads.hit5,
            ForwardAndReverseReads.hit6]

class StaggeredReads(object):
    """
    3 forward strand, 3 reverse strand reads.  Staggered starts and
    ends.

    FWD READS:                 gaa-t  a     (hit1)
                        a-ttacag att  aca   (hit2)
                       gatttaga             (hit3)
    REFERENCE:  (fwd)  GA TTACAG ATT  ACA
                 POS   01 234567 890  123
                (rev)  CT AATGTC TAA  TGT
    REV READS:          t aatctctt          (hit4)
                            tctctt-attt     (hit5)
                                   a ttg    (hit6)

    Note:
        - plurality insertion of A at 8;
        - plurality substitution C->G at 5;
        - plurality deletion of T at 9.
    """
    reference                  = "GATTACAGATTACA"
    expectedPluralityConsensus = "GATTAGAGAATACA"
    referenceWindow = (1, 0, 100)

    hit1 = AlignmentHitStub(7, FORWARD,
                            "G-ATTA",
                            "GAA-TA",
              InsertionQV=_("   -  "),
           SubstitutionQV=_("   -  "),
               DeletionQV=_("   -  "),
              DeletionTag=_("NNN-NN"),
                  MergeQV=_("   -  "))

    hit2 = AlignmentHitStub(1, FORWARD,
                            "ATTACAGATTACA",
                            "ATTACAGATTACA",
              InsertionQV=_("             "),
           SubstitutionQV=_("             "),
               DeletionQV=_("             "),
              DeletionTag=_("NNNNNNNNNNNNN"),
                  MergeQV=_("             "))

    hit3 = AlignmentHitStub(0, FORWARD,
                            "GA-TTACA",
                            "GATTTAGA",
              InsertionQV=_("  *     "),
           SubstitutionQV=_("        "),
               DeletionQV=_("        "),
              DeletionTag=_("NNNNNNNN"),
                  MergeQV=_("        "))

    hit4 = AlignmentHitStub(1, REVERSE,
                            "TAATGTC-T"[::-1],
                            "TAATCTCTT"[::-1],
              InsertionQV=_("         ")[::-1],
           SubstitutionQV=_("         ")[::-1],
               DeletionQV=_("         ")[::-1],
              DeletionTag=_("NNNNNNNNN")[::-1],
                  MergeQV=_("         ")[::-1])

    hit5 = AlignmentHitStub(4, REVERSE,
                            "TGTC-TAA--T"[::-1],
                            "TCTCTT-ATTT"[::-1],
              InsertionQV=_("      -    ")[::-1],
           SubstitutionQV=_("      -    ")[::-1],
               DeletionQV=_("      -    ")[::-1],
              DeletionTag=_("NNNNNN-NNNN")[::-1],
                  MergeQV=_("      -    ")[::-1])

    hit6 = AlignmentHitStub(10, REVERSE,
                            "A-TG"[::-1],
                            "ATTG"[::-1],
              InsertionQV=_("    ")[::-1],
           SubstitutionQV=_("    ")[::-1],
               DeletionQV=_("    ")[::-1],
              DeletionTag=_("NNNN")[::-1],
                  MergeQV=_("    ")[::-1])

    hits = [hit1, hit2, hit3, hit4, hit5, hit6]


class BigReads(object):
    """
    Large-ish hits, useful for stress-testing memory consumption.
    """
    length = 4000
    referenceWindow = (1, 0, length)
    seq = list("ACCT") * 1000
    np.random.seed(42)
    np.random.shuffle(seq)
    seq = "".join(seq)
    reference = expectedPluralityConsensus = seq

    middlingQV = np.array([7] * length, dtype=np.float32)
    hit1 = AlignmentHitStub(0, FORWARD, seq, seq,
                            InsertionQV    = middlingQV,
                            SubstitutionQV = middlingQV,
                            DeletionQV     = middlingQV,
                            DeletionTag    = _("N" * length),
                            MergeQV        = middlingQV)

    hits = [ hit1 ]

    @classmethod
    def manyHits(k, n):
        return [ AlignmentHitStub(0, FORWARD, k.seq, k.seq,
                                  InsertionQV    = k.middlingQV,
                                  SubstitutionQV = k.middlingQV,
                                  DeletionQV     = k.middlingQV,
                                  DeletionTag    = _("N" * k.length),
                                  MergeQV        = k.middlingQV)
                 for i in xrange(n) ]
