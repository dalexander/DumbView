
import numpy as np
from nose.tools import assert_equals


from DumbView.utils import Window, kSpannedIntervals, abut, holes


def test_intervals_1():
    """
    Intervals all covering the window
    """
    refWindow = Window(0, 100, 1010)
    start = np.array(np.array([100]*10, dtype=int), dtype=int)
    end   = np.array(np.array([110]*10, dtype=int), dtype=int)
    assert_equals([(100, 110)],
                   kSpannedIntervals(refWindow, 3, start, end))

def test_intervals_2():
    """
    Intervals not touching the window
    """
    refWindow = Window(0, 1, 10)
    start = np.array([0]*5 + [10]*5, dtype=int)
    end   = np.array([1]*5 + [15]*5, dtype=int)
    assert_equals([],
                  kSpannedIntervals(refWindow, 3, start, end))

def test_intervals_3():
    """
    Intervals covering the middle of the window -- "dromedary"
    """
    refWindow = Window(0, 0, 10)
    start = np.array([3]*10, dtype=int)
    end  = np.array([7]*10, dtype=int)
    assert_equals([(3, 7)],
                  kSpannedIntervals(refWindow, 3, start, end))

def test_intervals_4():
    """
    Two intervals at the fringes, with a hole in the middle --- "camel"
    """
    refWindow = Window(0, 100, 110)
    start = np.array([103]*5 + [107]*5, dtype=int)
    end   = np.array([105]*5 + [109]*5, dtype=int)
    assert_equals([(103,105), (107,109)],
                  kSpannedIntervals(refWindow, 3, start, end))


def test_intervals_5():
    """
    A case where there is nowhere 3-spanning coverage
    """
    refWindow = Window(0, 0, 10)
    reads = [ (x, x+1) for  x in xrange(0, 10) ]
    reads.append((0, 10))
    start, end = map(np.array, zip(*reads))
    assert_equals([ (x, x+1) for  x in xrange(0, 10) ],
                  kSpannedIntervals(refWindow, 2, start, end))
    assert_equals([],
                  kSpannedIntervals(refWindow, 3, start, end))


def test_abut():
    """
    Test abutting adjacent intervals
    """
    ints = [(s, s+1) for s in range(10)] + [(s, s+1) for s in range(20,30)]
    assert_equals([(0, 10), (20, 30)], abut(ints))


def test_holes_1():
    assert_equals([(0, 100)], holes(Window(0, 0, 100), []))

def test_holes_2():
    assert_equals([], holes(Window(0, 0, 100), [(0, 100)]))

def test_holes_3():
    """
    Holes for the dromedary test case
    """
    assert_equals([(0,3), (7,10)],
                  holes(Window(0, 0, 10), [(3,7)]))

def test_holes_4():
    """
    Holes for the camel test case
    """
    assert_equals([(3, 7)],
                  holes(Window(0, 0, 10), [(0,3), (7,10)]))
