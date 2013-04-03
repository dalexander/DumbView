
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
        start = refPos - 20
        end   = refPos + 20
    else:
        start, end = map(int, rest.split("-"))
    assert end >= start
    return Window(grok(key), start, end)

def windowFromGffString(windowString):
    w = windowFromString(windowString)
    return w._replace(start=w.start-1)

def subWindow(refWindow, subinterval):
    winId, winStart, winEnd = refWindow
    intS, intE = subinterval
    assert intS >= winStart
    assert intE <= winEnd
    return refWindow._replace(start=intS, end=intE)


# -------- Utility functions ---------

def integerValue(s):
    try:
        return int(s)
    except:
        return None

def grok(s):
    return integerValue(s) or s
