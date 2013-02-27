import numpy as np

ANSI_RED   = "\x1b[31m"
ANSI_RESET = "\x1b[39m"

def formatRed(useColor, txt):
    if useColor: return ANSI_RED + txt + ANSI_RESET
    else:        return txt

def formatReferenceCoordinates(refWindow):
    canvas = np.zeros(refWindow.end - refWindow.start, dtype="S1")
    canvas[:] = " "
    # Add one to convert to GFF coordinates...
    refCoords = np.arange(refWindow.start, refWindow.end) + 1
    for refPos in refCoords:
        if refPos % 20 == 0:
            sPos = str(refPos)
            sPosLen = len(sPos)
            cursorPos = refPos - sPosLen - refWindow.start
            if cursorPos > 0:
                canvas[cursorPos:(cursorPos+sPosLen)] = \
                    np.fromstring(sPos, dtype="S1")
    return canvas.tostring()

def formatSeparatorLine(refWindow):
    # Formats line like "--|----+----|----+--"
    canvas = np.zeros(refWindow.end - refWindow.start, dtype="S1")
    canvasCoords = np.arange(refWindow.start, refWindow.end) + 1
    canvas[:] = "="
    canvas[canvasCoords % 10 == 0] = "|"
    canvas[canvasCoords % 10 == 5] = "+"
    return canvas.tostring()

def formatAlignedRead(refWindow, cmpH5, rowNumber):
    canvas = np.zeros(refWindow.end - refWindow.start, dtype="S1")
    canvas[:] = " "
    clippedRead = cmpH5[rowNumber].clippedTo(refWindow.start, refWindow.end)
    fullRead = np.fromstring(clippedRead.read(orientation="genomic"), dtype="S1")
    transcript = np.fromstring(clippedRead.transcript(orientation="genomic"), dtype="S1")
    noInsertionsRead = np.extract(transcript != "I", fullRead)
    canvas[(clippedRead.tStart-refWindow.start):(clippedRead.tEnd-refWindow.start)] = \
        noInsertionsRead
    return canvas.tostring()

def formatUnalignedRead(refWindow, cmpH5, rowNumber, useColor=False):
    clippedRead = cmpH5[rowNumber].clippedTo(refWindow.start, refWindow.end)
    alnRead = clippedRead.read(orientation="genomic")
    transcript = clippedRead.transcript(orientation="genomic")
    output = ""
    for (readChar, transcriptChar) in zip(alnRead, transcript):
        if transcriptChar in "MR":  output += readChar
        elif transcriptChar in "I": output += formatRed(useColor, readChar.lower())
    return output

def formatAlignedReads(refWindow, cmpH5, rowNumbers):
    return [ formatAlignedRead(refWindow, cmpH5, rowNumber)
             for rowNumber in rowNumbers ]

def formatUnalignedReads(refWindow, cmpH5, rowNumbers, useColor=False):
    return [ formatUnalignedRead(refWindow, cmpH5, rowNumber, useColor)
             for rowNumber in rowNumbers ]

def formatWindow(cmpH5, refWindow, rowNumbers,
                 referenceTable=None, aligned=True, useColor=True):
    if referenceTable:
        referenceContig = referenceTable.byKey(cmpH5.referenceInfo(refWindow.name),
                                               refWindow.contigKey)
        referenceInWindow = referenceContig[refWindow.start:refWindow.end]
    else:
        referenceInWindow = None

    preMargin = " " * 10
    print preMargin + formatReferenceCoordinates(refWindow)
    if referenceInWindow:
        print "     Ref  " + referenceInWindow
    print preMargin + formatSeparatorLine(refWindow)

    if aligned:
        formattedReads = formatAlignedReads(refWindow, cmpH5, rowNumbers)
    else:
        formattedReads = formatUnalignedReads(refWindow, cmpH5, rowNumbers, useColor)

    for rn, ar in zip(rowNumbers, formattedReads):
        print ("%8d  " % rn)  + ar
