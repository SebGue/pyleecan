# -*- coding: utf-8 -*-

import numpy as np


def init_point(self, size_x, size_y, x, y):
    """
    Initialize the grid of point, of size size_x * size_y.
    And coordinates x,y (vectors).

    Parameters
    ----------
    size_x : integer
        Size of x.
    size_y : integer
        Size of y.
    x : nd-array, size: size_x (float)
        x coordinate.
    y : nd-array, size: size_x (float)
        y coordinate.

    Returns
    -------
    list_coord : nd-array, size: size_x*size_y x 2 (float)
        list of coordinate.

    """
    list_coord = np.zeros((size_x * size_y, 2))
    for i in range(x.size):
        for j in range(y.size):

            list_coord[x.size * j + i, :] = np.array([x[i], y[j]])

    return list_coord
