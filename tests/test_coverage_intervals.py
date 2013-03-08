
import numpy as np
from nose.tools import assert_equals

from DumbView.utils import find_k_spanned_intervals, abut


class Window(object):
    def __init__(self, start, end):
        self.start = start
        self.end  = end


def test_intervals_1():
    """
    Intervals all covering the window
    """
    refWindow = Window(100, 1010)
    start = np.array(np.array([100]*10, dtype=int), dtype=int)
    end   = np.array(np.array([110]*10, dtype=int), dtype=int)
    assert_equals([(100, 110)],
                   find_k_spanned_intervals(refWindow, 3, start, end))

def test_intervals_2():
    """
    Intervals not touching the window
    """
    refWindow = Window(1, 10)
    start = np.array([0]*5 + [10]*5, dtype=int)
    end   = np.array([1]*5 + [15]*5, dtype=int)
    assert_equals([],
                  find_k_spanned_intervals(refWindow, 3, start, end))

def test_intervals_3():
    """
    Intervals covering the middle of the window -- the "dromedary"
    test case
    """
    refWindow = Window(0, 10)
    start = np.array([3]*10, dtype=int)
    end  = np.array([7]*10, dtype=int)
    assert_equals([(3, 7)],
                  find_k_spanned_intervals(refWindow, 3, start, end))

def test_intervals_4():
    """
    Two intervals at the fringes, with a hole in the middle --- the
    "camel" test case
    """
    refWindow = Window(100, 110)
    start = np.array([103]*5 + [107]*5, dtype=int)
    end   = np.array([105]*5 + [109]*5, dtype=int)
    assert_equals([(103,105), (107,109)],
                  find_k_spanned_intervals(refWindow, 3, start, end))


def test_intervals_5():
    """
    A case where there is nowhere 3-spanning coverage
    """
    refWindow = Window(0, 10)
    reads = [ (x, x+1) for  x in xrange(0, 10) ]
    reads.append((0, 10))
    start, end = map(np.array, zip(*reads))
    assert_equals([ (x, x+1) for  x in xrange(0, 10) ],
                  find_k_spanned_intervals(refWindow, 2, start, end))
    assert_equals([],
                  find_k_spanned_intervals(refWindow, 3, start, end))



def test_abut():
    """
    Test abutting adjacent intervals
    """
    ints = [(s, s+1) for s in range(10)] + [(s, s+1) for s in range(20,30)]
    assert_equals([(0, 10), (20, 30)], abut(ints))
