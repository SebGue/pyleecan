# -*- coding: utf-8 -*-

import numpy as np


def init_mesh_BC(self, size_x, size_y, BC):
    """
    Init data structure to compute boundary condition in the sytem.

    Parameters
    ----------
    size_x : integer
        Number of horizontal points.
    size_y : integer
        Number of vertical point.
    BC : list of string (condition)
        List of boundary conditions.

    Returns
    -------
    BC_list : nd-array, size: m (int)
        Data strucuture.
    Periodic_point : nd-array, size: ? (int)
        Periodic point connection.

    """
    # BC: Boundary Condition, ["condi0","cond1","cond2","cond3"]
    # -> "HD": Homogenious Dicrichlet conditions (=0)
    # -> "P" : Periodic contitions
    # -> "AP": Anti Periodic Condition
    # BC_list -> Boudary condition

    nn2 = size_x * size_y
    ###########################################

    BC_list = np.zeros(nn2, dtype=np.uint8)
    Periodic_point = np.zeros(0, dtype=np.uint16)

    if np.all(BC == ["P", "HD", "P", "HD"]):

        # Initialyze the BC list

        # Vertical BC
        BC_list[::size_x] = 2
        BC_list[size_x - 1 :: size_x] = 2

        # Horizontal BC
        BC_list[:size_x] = 1
        BC_list[size_x * (size_y - 1) : nn2] = 1

        # Connexion of periodic points
        Periodic_point = -np.ones((size_y - 2, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(0, size_x * (size_y - 2), size_x) + size_x
        Periodic_point[:, 1] = Periodic_point[:, 0] + size_x - 1

    elif np.all(BC == ["AP", "HD", "AP", "HD"]):
        # Initialyze the BC list
        # Vertical BC
        BC_list[::size_x] = 2
        BC_list[size_x - 1 :: size_x] = 3
        # Horizontal BC
        BC_list[:size_x] = 1
        BC_list[size_x * (size_y - 1) : nn2] = 1

        Periodic_point = -np.ones((size_y - 2, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(0, size_x * (size_y - 2), size_x) + size_x
        Periodic_point[:, 1] = Periodic_point[:, 0] + size_x - 1

    elif np.all(BC == ["AP", "NA", "AP", "HD"]):
        # Initialyze the BC list
        # Vertical BC
        BC_list[::size_x] = 2
        BC_list[size_x - 1 :: size_x] = 3
        # Horizontal BC
        BC_list[size_x * (size_y - 1) : nn2] = 1

        Periodic_point = -np.ones((size_y - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(0, size_x * (size_y - 1), size_x)
        Periodic_point[:, 1] = Periodic_point[:, 0] + size_x - 1

    elif np.all(BC == ["AP", "HD", "AP", "NA"]):
        # Initialyze the BC list
        # Vertical BC
        BC_list[::size_x] = 2
        BC_list[size_x - 1 :: size_x] = 3

        # Horizontal BC
        BC_list[:size_x] = 1

        Periodic_point = -np.ones((size_y - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(size_x, size_x * (size_y), size_x)
        Periodic_point[:, 1] = Periodic_point[:, 0] + size_x - 1

    elif np.all(BC == ["P", "NA", "P", "HD"]):
        # Initialyze the BC list
        # Vertical BC
        BC_list[::size_x] = 2
        BC_list[size_x - 1 :: size_x] = 2
        # Horizontal BC
        BC_list[size_x * (size_y - 1) : nn2] = 1

        Periodic_point = -np.ones((size_y - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(0, size_x * (size_y - 1), size_x)
        Periodic_point[:, 1] = Periodic_point[:, 0] + size_x - 1

    elif np.all(BC == ["P", "HD", "P", "NA"]):
        # Initialyze the BC list
        # Vertical BC
        BC_list[::size_x] = 2
        BC_list[size_x - 1 :: size_x] = 2
        # Horizontal BC
        BC_list[:size_x] = 1

        Periodic_point = -np.ones((size_y - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(size_x, size_x * (size_y), size_x)
        Periodic_point[:, 1] = Periodic_point[:, 0] + size_x - 1

    else:
        raise NameError(
            "This boundary condition doesn't exist in this code or it's wrong"
        )

    return BC_list, Periodic_point
