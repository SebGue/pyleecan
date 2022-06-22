# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:08:06 2022

@author: Antoine Andrzejewski
"""
import numpy as np


def init_point(size_x, size_y, x, y):
    """
    Initaialyze grid of point, of size size_x * size_y.
    And coordinate x,y (vectors).

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


def init_permeabilty_cell(size_x, size_y, permeability_materials, cells_materials):
    """
    Initaialyze the permeabilty in each cell.

    Parameters
    ----------
    size_x : integer
        Number of horizontal points.
    size_y : integer
        Number of vertical point.
    permeability_materials : nd-array, size: (foat)
        permabilty of each materials.
    cells_materials : nd-array, size: m (int)
        materials in each cell.

    Returns
    -------
    permeability_cell : nd-array, size: n (int)
        Permeability in each cell.

    """
    permeability_cell = np.zeros((size_y - 1) * (size_x - 1), dtype=np.float64)
    permeability_cell[...] = permeability_materials[cells_materials - 1]

    return permeability_cell


def init_cell(size_x, size_y):
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


def init_mesh_BC(size_x, size_y, BC):
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

    BC_list = np.zeros(nn2, dtype=np.uint16)
    Periodic_point = np.zeros(0, dtype=np.uint16)

    if np.all(BC == ["HD", "HD", "HD", "HD"]):
        # Initialyze the BC list
        BC_list[1 : size_x - 1] = 1
        BC_list[size_x * (size_y - 1) + 1 : nn2 - 1] = 1
        BC_list[::size_y] = 1
        BC_list[size_x - 1 :: size_y] = 1

    elif np.all(BC == ["P", "HD", "P", "HD"]):
        # Initialyze the BC list
        # Vertical BC
        BC_list[::size_x] = 2
        BC_list[size_x - 1 :: size_x] = 2
        # Horizontal BC
        BC_list[:size_x] = 1
        BC_list[size_x * (size_y - 1) : nn2] = 1

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

        Periodic_point = -np.ones((size_y, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(0, size_x * (size_y), size_x)
        Periodic_point[:, 1] = Periodic_point[:, 0] + size_x - 1

    elif np.all(BC == ["P", "P", "P", "P"]):
        # Initialyze the BC list
        BC_list[1 : size_x - 1] = 2
        BC_list[size_x * (size_y - 1) + 1 : nn2 - 1] = 2

        BC_list[::size_y] = 2
        BC_list[size_x - 1 :: size_y] = 2

    elif np.all(BC == ["P", "AP", "P", "AP"]):
        # Initialyze the BC list
        # Vertical
        BC_list[::size_y] = 2
        BC_list[size_x - 1 :: size_y] = 3

        # Horizontal
        BC_list[:size_x] = 2
        BC_list[size_x * (size_y - 1) :] = 2
    else:
        print("wrong boundary conditions")

    return BC_list, Periodic_point


def numeroting_unknows(list_elem, BC_list, Periodic_point):
    """
    Pre-processing, numeroting the unknowns of the model in the linear system.

    Parameters
    ----------
    list_elem : nd-array, size: m (integers)
        tab of elements
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.
    Periodic_point : nd-array, size: ? (int)
        Periodic point connection.

    Returns
    -------
    Num_Unknowns : nd-array, size: n (int)
        Numerotation in the linear system.

    """
    Num_Unknowns = -np.ones(BC_list.size, dtype=np.int32)
    mask = BC_list == 0
    mask[Periodic_point[:, 0]] = True

    Num_Unknowns[mask] = np.arange(mask.sum(), dtype=np.int32)
    Num_Unknowns[Periodic_point[:, 1]] = Num_Unknowns[Periodic_point[:, 0]]

    return Num_Unknowns


def save_mesh(permeabiltiy_materials, Num_Unknowns, list_elem, x, y, BC_list):
    """
    SAve mesh

    Parameters
    ----------
    permeability_materials : nd-array, size: (foat)
        permabilty of each materials.
    Num_Unknowns : nd-array, size: n (int)
        Numerotation in the linear system.
    list_elem : nd-array, size: m (integers)
        tab of elements
    x : nd-array, size: size_x (float)
        x coordinate.
    y : nd-array, size: size_x (float)
        y coordinate.
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.

    Returns
    -------
    None.

    """
    f_handle = open("mesh.txt", "w")
    np.savetxt(
        f_handle,
        np.column_stack((list_elem, permeabiltiy_materials)),
        fmt="%u",
        header=str(permeabiltiy_materials.size) + " " + str(Num_Unknowns.size),
    )
    list_coord = np.zeros((Num_Unknowns.size, 2))
    for i in range(x.size):
        for j in range(y.size):
            list_coord[x.size * j + i, :] = np.array([x[i], y[j]])
    f_handle.close()
    f_handle = open("mesh.txt", "a")
    np.savetxt(
        f_handle,
        np.column_stack((list_coord, Num_Unknowns, BC_list)),
        fmt="%f %f %d %d",
    )
    f_handle.close()
    print("mesh is saved")


def init_reluc(list_elem, list_coord, mu0, la, mode):
    """
    Parameters
    ----------
    list_elem : nd-array of integers ( m x 4 )
        list of the element
    list_coord : nd-array of coordinate ( n x 2 )
        coordinate of each vertices

    Returns
    -------
    Reluc_list : Reluc_list (m x 4 )
        Reluctance of each elements
    """
    R_reluc = np.zeros((list_elem.shape[0], 4))
    if mode == "cartesian":
        # length and width
        h_x = np.linalg.norm(
            list_coord[list_elem[:, 0]] - list_coord[list_elem[:, 1]], axis=1, ord=2
        )
        h_y = np.linalg.norm(
            list_coord[list_elem[:, 1]] - list_coord[list_elem[:, 2]], axis=1, ord=2
        )

        R_reluc[:, 0] = 0.5 / (mu0 * la) * h_y / h_x
        R_reluc[:, 1] = 0.5 / (mu0 * la) * h_x / h_y
        R_reluc[:, 2] = 0.5 / (mu0 * la) * h_y / h_x
        R_reluc[:, 3] = 0.5 / (mu0 * la) * h_x / h_y

    elif mode == "polar":
        theta = np.abs(list_coord[list_elem[:, 0], 0] - list_coord[list_elem[:, 1], 0])
        R0 = list_coord[list_elem[:, 0], 1]
        R1 = list_coord[list_elem[:, 3], 1]

        R_MEAN = (R0 + R1) / 2

        R_reluc[:, 0] = (np.log(R1 / R_MEAN)) / (mu0 * theta * la)
        R_reluc[:, 1] = 0.5 * theta / (mu0 * la * np.log(R1 / R0))
        R_reluc[:, 2] = (np.log(R_MEAN / R0)) / (mu0 * theta * la)
        R_reluc[:, 3] = 0.5 * theta / (mu0 * la * np.log(R1 / R0))
    else:
        print("Wrong mode, choice between ")

    return R_reluc
