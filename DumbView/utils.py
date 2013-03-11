
import numpy as np
from pbcore.io import ReferenceTable
from pbcore.io.rangeQueries import projectIntoRange

class SuperReferenceTable(ReferenceTable):
    def byKey(self, key):
        """
        Lookup key by: (1) name, (2) localId (3) localId as int if int-convertible
        """
        try: return self.byName(key)
        except: pass
        try: return self.byLocalId(key)
        except: pass
        try: return self.byLocalId(int(key))
        except:
            raise KeyError("Contig not found")

    def keyToLocalId(self, key):
        return self.byKey(key).localId


class Window(object):
    def __init__(self, contigKey, start, end):
        self.contigKey  = contigKey
        self.start = start
        self.end   = end

    def __repr__(self):
        return "%s:%d-%d" % (self.contigKey, self.start, self.end)

    def __len__(self):
        return (self.end - self.start)

    def __eq__(self, other):
        return (self.contigKey == other.contigKey and
                self.start     == other.start and
                self.end       == other.end)

    @staticmethod
    def fromString(windowString):
        splitOnColon = windowString.split(":")
        if len(splitOnColon) < 2:
            key = 1
        else:
            key, rest = splitOnColon
        if "-" not in rest:
            refPos = int(rest)
            start = refPos - 10
            end   = refPos + 10
        else:
            start, end = map(int, rest.split("-"))
        assert end >= start
        return Window(key, start, end)

    @staticmethod
    def fromGffString(windowString):
        w = Window.fromString(windowString)
        w.start -= 1
        return w

    def subWindow(self, start, end):
        assert start >= self.start
        assert end   <= self.end
        return Window(self.contigKey, start, end)

def readsInWindow(cmpH5, window, depthLimit=None, minMapQV=0, strategy="fileorder"):
    """
    Return up to `depthLimit` reads (as row numbers integers) where
    the mapped reference intersects the window.  If depthLimit is None,
    return all the reads meeting the criteria.

    `strategy` can be:
      - "longest" --- get the reads with the longest length in the window
      - "spanning" --- get only the reads spanning the window
      - "fileorder" --- get the reads in file order

    """
    assert strategy in {"longest", "spanning", "fileorder"}

    def depthCap(lst):
        if depthLimit is not None: return lst[:depthLimit]
        else: return lst

    contigLocalId = cmpH5.referenceInfo(window.contigKey).ID

    rowNumbers = cmpH5.readsInRange(contigLocalId, window.start, window.end,
                                    justIndices=True)

    rowNumbers = rowNumbers[cmpH5.MapQV[rowNumbers] >= minMapQV]

    if strategy == "fileorder":
        return depthCap(rowNumbers)

    tStartTruncated = np.maximum(window.start, cmpH5.tStart[rowNumbers])
    tEndTruncated   = np.minimum(window.end,   cmpH5.tEnd[rowNumbers])
    lengthsInWindow = tEndTruncated - tStartTruncated

    if strategy == "spanning":
        return depthCap(rowNumbers[lengthsInWindow==(window.end-window.start)])
    elif strategy == "longest":
        ordering = np.lexsort((rowNumbers, -lengthsInWindow))
        return depthCap(rowNumbers[ordering])



def kSpannedIntervals(refWindow, k, start, end):
    """
    Find intervals in the window that are k-spanned by the reads.

    Given:
     `refWindow`: the window under consideration
     `k`: the number of reads that must span intervals to be returned
     `start`, `end`: numpy arrays of start and end coordinates for reads,
       where the extent of each read is [start, end).  Must be ordered
       so that `start` is sorted in ascending order.

    Find a maximal set of maximal disjoint intervals within
    refWindow such that each interval is spanned by at least k reads.
    Intervals are returned in sorted order, as a list of (start, end)
    tuples.

    Note that this is a greedy search procedure and may not always
    return the optimal solution, in some sense.  However it will
    always return the optimal solutions in the most common cases.
    """
    assert k >= 1


    # Translate the start, end to coordinate system where
    # refWindow.start is 0.
    start = start - refWindow.start
    end   = end - refWindow.start
    winStart = 0
    winEnd   = refWindow.end - refWindow.start

    # Truncate to bounds implied by refWindow
    start = np.maximum(winStart, start)
    end   = np.minimum(winEnd,   end)

    positions = np.arange(winEnd - winStart, dtype=int)
    coverage = projectIntoRange(start, end,
                                winStart, winEnd)
    x = -1
    y = 0
    intervalsFound = []

    while y < winEnd:
        # Step 1: let x be the first pos >= y that is k-covered
        eligible = np.flatnonzero((positions >= y) & (coverage >= k))
        if len(eligible) > 0:
            x = eligible[0]
        else:
            break

        # Step 2: extend the window [x, y) until [x, y) is no longer
        # k-spanned.  Do this by setting y to the k-th largest `end`
        # among reads covering x
        eligible = end[(start <= x)]
        eligible.sort()
        if len(eligible) >= k:
            y = eligible[-k]
        else:
            break

        intervalsFound.append((x, y))

    # Translate intervals back
    return [ (s + refWindow.start,
              e + refWindow.start) for (s, e) in intervalsFound ]


def abut(intervals):
    """
    Abut adjacent intervals.  Useful for debugging...
    """
    output = []
    lastS = None
    lastE = None
    for (s, e) in intervals:
        if s == lastE:
            lastS, lastE = lastS, e
        else:
            if lastS is not None:
                output.append((lastS, lastE))
            lastS, lastE = s, e
    output.append((lastS, lastE))
    return output


def holes(refWindow, intervals):
    """
    Given a window and a set of disjoint subintervals, return the
    "holes", which are the intervals of the refWindow not covered by
    the given subintervals.
    """
    output = []
    intervals = sorted(intervals)
    lastE = refWindow.start
    for (s, e) in intervals:
        if s > lastE:
            output.append((lastE, s))
        lastE = e
    if lastE < refWindow.end:
        output.append((lastE, refWindow.end))
    return output
