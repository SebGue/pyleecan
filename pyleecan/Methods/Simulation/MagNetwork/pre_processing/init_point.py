# -*- coding: utf-8 -*-

import numpy as np


def init_point(self, N_point_theta, N_point_r, theta, r):
    """
    Initialize the grid of points (size : N_point_theta * N_point_r) in the r and theta directions

    Parameters
    ----------
    N_point_theta : integer
        Discretization point in the theta direction
    N_point_r : integer
        Discretization points in the r direction
    theta : nd-array (size: N_point_theta of floats)
        Coordinates of points in the theta direction
    r : nd-array (size: N_point_r of floats)
        Coordinates of points in the r direction

    Returns
    -------
    list_coord : nd-array (size: N_point_r * N_point_theta * 2 of loats )
        list of point coordinates

    """
    if N_point_theta != theta.size:
        raise NameError("Wrong number of points in theta direction, N_point_theta= "+str(N_point_theta)+" and theta.size= "+str(theta.size)+" must be equal")
    if N_point_r != r.size:
        raise NameError("Wrong number of points in radial direction.")

    list_coord = np.zeros((N_point_theta * N_point_r, 2))

    list_coord[:, 0] = np.tile(theta, N_point_r)
    list_coord[:, 1] = np.repeat(r, N_point_theta)

    return list_coord
