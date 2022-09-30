# -*- coding: utf-8 -*-

import numpy as np


def init_cell(self, N_point_theta, N_point_r):
    """
    Compute the list of the elements in regular grid.

    Parameters
    ----------
    N_point_theta : integer
        Number of horizontal points.
    N_point_r : integer
        Number of vertical point.

    Returns
    -------
    list_elem : nd-array, size: m (integers)
        tab of rectangular elements

    """

    list_elem = np.zeros(((N_point_theta - 1) * (N_point_r - 1), 4), dtype=np.uint16)
    for i in range(N_point_r - 1):
        for j in range(N_point_theta - 1):
            element = (N_point_theta - 1) * i + j
            list_elem[element, 0] = N_point_theta * i + j
            list_elem[element, 1] = N_point_theta * i + j + 1
            list_elem[element, 2] = N_point_theta * (i + 1) + j + 1
            list_elem[element, 3] = N_point_theta * (i + 1) + j
    return list_elem
