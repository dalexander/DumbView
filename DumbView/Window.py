
from collections import namedtuple

Window = namedtuple("Window", ["refId", "start", "end"])

def windowFromString(windowString):
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

def windowFromGffString(windowString):
    w = windowFromString(windowString)
    return w._replace(start=w.start-1)
