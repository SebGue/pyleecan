# -*- coding: utf-8 -*-

import numpy as np


def init_cell(self, size_x, size_y):
    """
    Compute the list of the elements in regular grid.

    Parameters
    ----------
    size_x : integer
        Number of horizontal points.
    size_y : integer
        Number of vertical point.

    Returns
    -------
    list_elem : nd-array, size: m (integers)
        tab of rectangular elements

    """

    list_elem = np.zeros(((size_x - 1) * (size_y - 1), 4), dtype=np.uint16)
    for i in range(size_y - 1):
        for j in range(size_x - 1):
            element = (size_x - 1) * i + j
            list_elem[element, 0] = size_x * i + j
            list_elem[element, 1] = size_x * i + j + 1
            list_elem[element, 2] = size_x * (i + 1) + j + 1
            list_elem[element, 3] = size_x * (i + 1) + j
    return list_elem
