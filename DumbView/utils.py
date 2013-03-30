
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

    localRefId = cmpH5.referenceInfo(window.refId).ID
    rowNumbers = cmpH5.readsInRange(localRefId, window.start, window.end,
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
