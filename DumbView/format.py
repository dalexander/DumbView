# -*- coding: utf-8 -*-

SPARKS = u' ▁▂▃▄▅▆▇'

import numpy as np

from itertools import groupby
from operator import itemgetter

from .consensus import consensus, align

ANSI_RED     = "\x1b[31m"
ANSI_GREEN   = "\x1b[32m"
ANSI_REVERSE = "\x1b[7m"
ANSI_RESET   = "\x1b[0m"

def red(txt):
    return ANSI_RED + txt + ANSI_RESET

def green(txt):
    return ANSI_GREEN + txt + ANSI_RESET

def reverseRed(txt):
    return ANSI_RED + ANSI_REVERSE + txt + ANSI_RESET

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

def formatAlignedRead(cmpH5, refWindow, rowNumber):
    canvas = np.zeros(refWindow.end - refWindow.start, dtype="S1")
    canvas[:] = " "
    clippedRead = cmpH5[rowNumber].clippedTo(refWindow.start, refWindow.end)
    fullRead = np.fromstring(clippedRead.read(orientation="genomic"), dtype="S1")
    transcript = np.fromstring(clippedRead.transcript(orientation="genomic"), dtype="S1")
    noInsertionsRead = np.extract(transcript != "I", fullRead)
    canvas[(clippedRead.tStart-refWindow.start):(clippedRead.tEnd-refWindow.start)] = \
        noInsertionsRead
    return canvas.tostring()

def formatAlignedRead3(cmpH5, refWindow, referenceInWindow, rowNumber, useColor=False):
    try:
        clippedRead = cmpH5[rowNumber].clippedTo(refWindow.start, refWindow.end)
    except:
        return ""

    read = clippedRead.read(orientation="genomic")
    transcript = clippedRead.transcript(orientation="genomic")
    validMoves = set(["R", "D", "M"])

    # get read absent any insertions (should be of refWindow length)
    read, transcript = [
        "".join(itr) for itr in zip(*[
            (r, x) for r, x in zip(read, transcript)
            if x in validMoves])]

    assert len(read) == len(referenceInWindow)

    # if you have a sequence of Ds and Ms, move Ds up
    def bubbleUpGaps(transcript):
        validMoves = set(["D", "M"])
        state = []
        for x in transcript:
            if x in validMoves:
                state.append(x)
            else:
                state.sort()
                for s in state:
                    yield s
                state = []
                yield x
        state.sort()
        for s in state:
            yield s

    # group by homopolymer in the reference, bubble up gaps to head of the HP
    transcript = "".join(
        "".join(bubbleUpGaps(x for _, x in rx))
        for _, rx in groupby(zip(referenceInWindow, transcript), itemgetter(0)))

    startGap = clippedRead.tStart - refWindow.start

    def renderRead():
        for _ in xrange(startGap):
            yield " "

        for x, r in zip(transcript, read):
            if x == "R":
                yield reverseRed(r) if useColor else r
            elif x == "D":
                yield "-"
            elif x == "M":
                yield r

    return "".join(renderRead())

def formatAlignedRead2(cmpH5, refWindow, rowNumber, useColor=False):
    try:
        clippedRead = cmpH5[rowNumber].clippedTo(refWindow.start, refWindow.end)
    except:
        return ""

    read = clippedRead.read(orientation="genomic")
    transcript = clippedRead.transcript(orientation="genomic")
    rendered = ""
    for x, r in zip(transcript, read):
        if x == "R":
            rendered += reverseRed(r) if useColor else r
        elif x == "D":
            rendered += "-"
        elif x == "M":
            rendered += r
    startGap = clippedRead.tStart - refWindow.start
    return " "*startGap + rendered


def formatUnalignedRead(cmpH5, refWindow, rowNumber, useColor=False):
    # FIXME!  This code is incorrect for reads that start partway through the window!
    # Ex, reads 14 and 1784 from ref000001:1-100 of job 038537
    clippedRead = cmpH5[rowNumber].clippedTo(refWindow.start, refWindow.end)
    alnRead = clippedRead.read(orientation="genomic")
    transcript = clippedRead.transcript(orientation="genomic")
    output = ""
    for (readChar, transcriptChar) in zip(alnRead, transcript):
        if transcriptChar in "MR":  output += readChar
        elif transcriptChar in "I":
            if useColor: output += red(readChar.lower())
            else:        output += readChar.lower()
    return output

def formatAlignedReads2(cmpH5, refWindow, referenceInWindow, rowNumbers, useColor=False):
    return [ formatAlignedRead3(cmpH5, refWindow, referenceInWindow, rowNumber, useColor)
             for rowNumber in rowNumbers ]

def formatAlignedReads(cmpH5, refWindow, rowNumbers, useColor=False):
    return [ formatAlignedRead2(cmpH5, refWindow, rowNumber, useColor)
             for rowNumber in rowNumbers ]

def formatUnalignedReads(cmpH5, refWindow, rowNumbers, useColor=False):
    return [ formatUnalignedRead(cmpH5, refWindow, rowNumber, useColor)
             for rowNumber in rowNumbers ]

def formatConsensus(cmpH5, refWindow, rowNumbers, refTable):
    cssObj = consensus(cmpH5, refWindow, rowNumbers, refTable)
    print "     CSS  " + cssObj.sequence
    print " " * 10 + spark(cssObj.confidence)

def formatWindow(cmpH5, refWindow, rowNumbers,
                 referenceTable=None, aligned=True, useColor=True, consensus=True):

    if referenceTable:
        refName = cmpH5.referenceInfo(refWindow.refId).FullName
        referenceInWindow = referenceTable[refName].sequence[refWindow.start:refWindow.end]
    else:
        referenceInWindow = None

    preMargin = " " * 10
    print preMargin + formatReferenceCoordinates(refWindow)
    if referenceInWindow:
        print "     Ref  " + referenceInWindow
    print preMargin + formatSeparatorLine(refWindow)

    if aligned and referenceTable:
        formattedReads = formatAlignedReads2(cmpH5, refWindow, referenceInWindow, rowNumbers, useColor)
    elif aligned:
        formattedReads = formatAlignedReads(cmpH5, refWindow, rowNumbers, useColor)
    else:
        formattedReads = formatUnalignedReads(cmpH5, refWindow, rowNumbers, useColor)

    for rn, ar in zip(rowNumbers, formattedReads):
        print ("%8d  " % rn)  + ar

    if referenceTable and consensus:
        print
        print preMargin + formatReferenceCoordinates(refWindow)
        print preMargin + formatSeparatorLine(refWindow)
        formatReferenceAndConsensus(cmpH5, refWindow, referenceTable, rowNumbers)

def spark(arr):
    idx = (np.array(arr, dtype=np.uint)/8).clip(0, len(SPARKS)-1)
    #print idx
    return "".join(SPARKS[i] for i in idx)

def formatReferenceAndConsensus(cmpH5, refWindow, refTable, rowNumbers):
    refName = cmpH5.referenceInfo(refWindow.refId).FullName
    cssObj = consensus(cmpH5, refWindow, refTable, rowNumbers)
    refInWindow = refTable[refName].sequence[refWindow.start:refWindow.end]
    alnRef, transcript, alnQuery = align(refInWindow, cssObj.sequence)

    refChars = []
    cssChars = []

    for (r, t, q) in zip(alnRef, transcript, alnQuery):
        if t == "M":
            refChars.append(r)
            cssChars.append(q)
        elif t == "R":
            refChars.append(red(r))
            cssChars.append(green(q))
        elif t == "I":
            cssChars.append(green(q))
        elif t == "D":
            refChars.append(red(r))

    print "     ref  " + "".join(refChars)
    print "     css  " + "".join(cssChars)
    print " " * 10 + spark(cssObj.confidence)


def formatIndividualAlignments(cmpH5, refWindow, rowNumbers):
    for row in rowNumbers:
        print "--"
        print "Row number %d" % row
        print cmpH5[row].clippedTo(refWindow.start, refWindow.end)
