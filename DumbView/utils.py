
import numpy as np, sys
from pbcore.io import ReferenceTable
from pbcore.io.rangeQueries import projectIntoRange

def die(msg):
    print msg
    sys.exit(-1)

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
