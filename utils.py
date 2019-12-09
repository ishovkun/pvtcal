#!/usr/bin/env python
from enum import Enum
import numpy as np
log = np.log

def interp_lin(x1, x2, y1, y2, x):
    """
    simple linear interpolation
    y1 and y2 are the function values in the points
    x1 and x2, respectively
    x is the point where we want to evaluate the
    function
    """
    return (y2 - y1) / (x2 - x1) * (x - x1) + y1

def interp_log(x1, x2, y1, y2, x):
    """
    simple logarithmic interpolation
    no checks
    """
    return np.exp( interp_lin(x1, x2, log(y1), log(y2), x) )

class Position(Enum):
    """
    to check whether we are between
    undersaturated branches or out of region
    """
    BETWEEN = 1
    ABOVE = 2
    BELOW = 3


def binary_search(value, array, l, h):
    '''
    regular binary search
    '''
    if (h - l == 1):
        if (value == array[h]):
            return h
        elif (value == array[l]):
            return l
        else:
            raise KeyError("element not found")

    m = int((h - l) / 2)
    print(m, array[m])
    if (value > array[m]):
        return binary_search(value, array, m, h)
    elif (value < array[m]):
        return binary_search(value, array, l, m)
    else: #if ==
        return m

def binary_search_extr(value, array, l, h, shift):
    '''
    a version of binary search that provides two
    closest indices for interpolation/extrapolation
    Complexity: O(log(size(array)))
    '''
    if (h - l == 1):
        return l, h

    m = int(l + (h - l) / 2)
    print(l,m,h)
    if (value > array[m] - shift):
        return binary_search_extr(value, array, m, h, shift)
    elif (value < array[m]):
        return binary_search_extr(value, array, l, m, shift)
    else: #if ==
        if (array[h] - array[m] > array[m] - array[l]):
            return l, m
        else: return m, h

def binarySearchExtrapolationIndices(value, array, shift = 0):
    '''
    value is a single value (float)
    find surrounding undersaturated branches
    shift is to search relative pressure
    Complexity: O(log(size(array)))
    '''
    assert sorted(array), "Array is not sorted"
    return binary_search_extr(value, array, 0, len(array)-1, shift)

def findSurroundingElements(value, array, shift = 0):
    """
    value is a single value (float)
    find surrounding undersaturated branches
    shift is to search relative pressure
    """
    lower = 0
    for i in range(len(array)):
        if (array[i] - shift <= value):
            lower = i

    upper = len(array) - 1
    for i in range(len(array)-1, -1, -1):
        if (array[i] - shift > value):
            upper = i

    if (value <= array[0] - shift):
        return Position.BELOW, lower, upper
    elif (value > array[-1] - shift):
        return Position.ABOVE, lower, upper
    else:
        return Position.BETWEEN, lower, upper


def curves_intersect(x1, x2, y1, y2):
    """
    return true if two curves y1 = y1(x1) and y2 = y2(x2)
    intersect
    x1, x2, y1, and y2 are all either lists or np arrays
    """
    assert(len(x1) == len(y1))
    assert(len(x2) == len(y2))
    # split curves into segments : y = a1 x + b1 and y = a2 x + b2
    # and check all pairs for intersection
    for i in range(1, len(x1)):
        # build segment line y = a1 x + b1
        assert x1[i] != x1[i-1]
        a1 = (y1[i] - y1[i-1]) / (x1[i] - x1[i-1])
        b1 = y1[i] - a1 * x1[i]
        for j in range(1, len(x2)):
            # second line : y = a2 x + b2
            assert x2[j] != x2[j-1]
            a2 = (y2[j] - y2[j-1]) / (x2[j] - x2[j-1])
            b2 = y2[j] - a2 * x2[j]

            # compute intersection
            if (abs(a1 - a2) < 1e-6):
                if (abs (b1 - b2) < 1e-6):
                    return True # conicide
                else:
                    continue

            # intersection
            x = - (b2 - b1) / (a2 - a1)
            if ( x >= x1[i-1] and x <= x1[i] and x >= x2[j-1] and x <= x2[j]):
                return True
    return False
