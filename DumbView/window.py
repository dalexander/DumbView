
from collections import namedtuple

Window = namedtuple("Window", ["refId", "start", "end"])

def windowFromString(windowString):
    splitOnColon = windowString.split(":")
    if len(splitOnColon) < 2:
        key = 1
        rest = windowString
    else:
        key, rest = splitOnColon
    if "-" not in rest:
        refPos = int(rest)
        start = refPos
        end   = None
    else:
        start, end = map(int, rest.split("-"))
    assert end is None or end >= start
    return Window(grok(key), start, end)

def windowFromGffString(windowString):
    w = windowFromString(windowString)
    return w._replace(start=w.start-1)

def windowsFromGffStrings(gffStrings):
    return map(windowFromGffString, gffStrings.split(","))

def windowToGffString(window):
    return "%s:%d-%d" % (window.refId, window.start+1, window.end)

def subWindow(refWindow, subinterval):
    winId, winStart, winEnd = refWindow
    intS, intE = subinterval
    assert intS >= winStart
    assert intE <= winEnd
    return refWindow._replace(start=intS, end=intE)

def makeDisplayWindow(contigLen, width, refWindow):
    refId, refStart, refEnd = refWindow
    if refEnd == None:
        refEnd = refStart + width//2
        refStart = refStart - width//2
    refStart = max(0, refStart)
    refEnd   = min(contigLen, refEnd)
    return Window(refId, refStart, refEnd)


# -------- Utility functions ---------

def integerValue(s):
    try:
        if s.isdigit():
            return int(s)
    except:
        pass
    return None

def grok(s):
    iv = integerValue(s)
    if iv is None:
        return s
    else:
        return iv
