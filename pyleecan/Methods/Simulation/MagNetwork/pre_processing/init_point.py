# -*- coding: utf-8 -*-

import numpy as np


def init_point(self, N_point_theta, N_point_r, theta, r):
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
    # if N_point_theta != theta.size:
    #     raise NameError("Wrong number of points in theta direction.")
    if N_point_r != r.size:
        raise NameError("Wrong number of points in radial direction.")

    list_coord = np.zeros((N_point_theta * N_point_r, 2))
    # for i in range(N_point_theta):
    #     for j in range(N_point_r):
    #         list_coord[N_point_theta * j + i, :] = np.array([r[i], theta[j]])

    list_coord[:, 0] = np.tile(theta, N_point_r)
    list_coord[:, 1] = np.repeat(r, N_point_theta)

    return list_coord
