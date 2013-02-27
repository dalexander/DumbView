
import numpy as np
from pbcore.io import ReferenceTable

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

def readsInWindow(cmpH5, window, depth, minMapQV=0, strategy="longest"):
    """
    Return up to `depth` reads (as row numbers integers) where the mapped
    reference intersects the window.

    `strategy` can be:
      - "longest" --- get the reads with the longest length in the window
      - "spanning" --- get only the reads spanning the window
      - "fileorder" --- get the reads in file order

    """
    assert strategy in {"longest", "spanning", "fileorder"}
    contigLocalId = cmpH5.referenceInfo(window.contigKey).ID

    # FIXME: this is slowish ... need to use readLocator nonsense
    rowNumbers = cmpH5.readsInRange(contigLocalId, window.start, window.end,
                                    justIndices=True)

    rowNumbers = rowNumbers[cmpH5.MapQV[rowNumbers] > minMapQV]

    if strategy == "fileorder":
        return rowNumbers[:depth]

    tStartTruncated = np.maximum(window.start, cmpH5.tStart[rowNumbers])
    tEndTruncated   = np.minimum(window.end,   cmpH5.tEnd[rowNumbers])
    lengthsInWindow = tEndTruncated - tStartTruncated

    if strategy == "spanning":
        return rowNumbers[lengthsInWindow==(window.end-window.start)][:depth]
    elif strategy == "longest":
        ordering = np.lexsort((rowNumbers, -lengthsInWindow))
        return rowNumbers[ordering[:depth]]
