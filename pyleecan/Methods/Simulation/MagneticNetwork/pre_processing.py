# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:08:06 2022

@author: Antoine Andrzejewski
"""
import numpy as np


def init_point(self, size_x, size_y, x, y):
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
    list_coord = np.zeros((self.size_x * self.size_y, 2))
    for i in range(self.x.size):
        for j in range(self.y.size):
            list_coord[self.x.size * j + i, :] = np.array([self.x[i], self.y[j]])

    return list_coord


def init_permeabilty_cell(
    self, size_x, size_y, permeability_materials, cells_materials
):
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
    permeability_cell = np.zeros(
        (self.size_y - 1) * (self.size_x - 1), dtype=np.float64
    )
    permeability_cell[...] = self.permeability_materials[self.cells_materials - 1]

    return permeability_cell


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

    list_elem = np.zeros(((self.size_x - 1) * (self.size_y - 1), 4), dtype=np.uint16)
    for i in range(self.size_y - 1):
        for j in range(self.size_x - 1):
            element = (self.size_x - 1) * i + j
            list_elem[element, 0] = self.size_x * i + j
            list_elem[element, 1] = self.size_x * i + j + 1
            list_elem[element, 2] = self.size_x * (i + 1) + j + 1
            list_elem[element, 3] = self.size_x * (i + 1) + j
    return list_elem


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

    nn2 = self.size_x * self.size_y
    ###########################################

    BC_list = np.zeros(nn2, dtype=np.uint8)
    Periodic_point = np.zeros(0, dtype=np.uint16)

    if np.all(self.BC == ["HD", "HD", "HD", "HD"]):
        # Initialyze the BC list
        BC_list[1 : self.size_x - 1] = 1
        BC_list[self.size_x * (self.size_y - 1) + 1 : nn2 - 1] = 1
        BC_list[:: self.size_y] = 1
        BC_list[self.size_x - 1 :: self.size_y] = 1

        Periodic_point = -np.ones((size_y - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(
            self.size_x, self.size_x * (self.size_y), self.size_x
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + self.size_x - 1

    elif np.all(self.BC == ["P", "NA", "P", "HD"]):
        # Initialyze the BC list
        # Vertical BC
        BC_list[:: self.size_x] = 2
        BC_list[self.size_x - 1 :: self.size_x] = 2
        # Horizontal BC
        BC_list[self.size_x * (self.size_y - 1) : nn2] = 1

        Periodic_point = -np.ones((self.size_y - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(
            0, self.size_x * (self.size_y - 1), self.size_x
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + self.size_x - 1

    elif np.all(self.BC == ["P", "HD", "P", "NA"]):
        # Initialyze the BC list
        # Vertical BC
        BC_list[:: self.size_x] = 2
        BC_list[self.size_x - 1 :: self.size_x] = 2
        # Horizontal BC
        BC_list[: self.size_x] = 1

        Periodic_point = -np.ones((self.size_y - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(
            self.size_x, self.size_x * (self.size_y), self.size_x
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + self.size_x - 1

    elif np.all(self.BC == ["P", "P", "P", "P"]):
        # Initialyze the BC list
        BC_list[1 : self.size_x - 1] = 2
        BC_list[self.size_x * (self.size_y - 1) + 1 : nn2 - 1] = 2

        BC_list[:: self.size_y] = 2
        BC_list[self.size_x - 1 :: self.size_y] = 2

    elif np.all(self.BC == ["P", "AP", "P", "AP"]):
        # Initialyze the BC list
        # Vertical
        BC_list[:: self.size_y] = 2
        BC_list[self.size_x - 1 :: self.size_y] = 3

        # Horizontal
        BC_list[: self.size_x] = 2
        BC_list[self.size_x * (self.size_y - 1) :] = 2
    else:
        print("wrong boundary conditions")

    return BC_list, Periodic_point


def numeroting_unknows(self, list_elem, BC_list, Periodic_point):
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
    Num_Unknowns = -np.ones(self.BC_list.size, dtype=np.int32)
    mask = self.BC_list == 0
    mask[self.Periodic_point[:, 0]] = True

    Num_Unknowns[mask] = np.arange(mask.sum(), dtype=np.int32)
    Num_Unknowns[self.Periodic_point[:, 1]] = Num_Unknowns[self.Periodic_point[:, 0]]

    return Num_Unknowns


def save_mesh(self, permeabiltiy_materials, Num_Unknowns, list_elem, x, y, BC_list):
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
        np.column_stack((self.list_elem, self.permeabiltiy_materials)),
        fmt="%u",
        header=str(self.permeabiltiy_materials.size)
        + " "
        + str(self.Num_Unknowns.size),
    )
    list_coord = np.zeros((self.Num_Unknowns.size, 2))
    for i in range(self.x.size):
        for j in range(self.y.size):
            list_coord[self.x.size * j + i, :] = np.array([self.x[i], self.y[j]])
    f_handle.close()
    f_handle = open("mesh.txt", "a")
    np.savetxt(
        f_handle,
        np.column_stack((list_coord, self.Num_Unknowns, self.BC_list)),
        fmt="%f %f %d %d",
    )
    f_handle.close()
    print("mesh is saved")


def init_reluc(self, list_elem, list_coord, mu0, la, mode):
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
    R_reluc = np.zeros((self.list_elem.shape[0], 4))
    if self.mode == "cartesian":
        # length and width
        h_x = np.linalg.norm(
            self.list_coord[self.list_elem[:, 0]]
            - self.list_coord[self.list_elem[:, 1]],
            axis=1,
            ord=2,
        )
        h_y = np.linalg.norm(
            self.list_coord[self.list_elem[:, 1]]
            - self.list_coord[self.list_elem[:, 2]],
            axis=1,
            ord=2,
        )

        R_reluc[:, 0] = 0.5 / (self.mu0 * self.la) * h_y / h_x
        R_reluc[:, 1] = 0.5 / (self.mu0 * self.la) * h_x / h_y
        R_reluc[:, 2] = 0.5 / (self.mu0 * self.la) * h_y / h_x
        R_reluc[:, 3] = 0.5 / (self.mu0 * self.la) * h_x / h_y

    elif self.mode == "polar":
        theta = np.abs(
            self.list_coord[self.list_elem[:, 0], 0]
            - self.list_coord[self.list_elem[:, 1], 0]
        )
        R0 = self.list_coord[self.list_elem[:, 0], 1]
        R1 = self.list_coord[self.list_elem[:, 3], 1]

        R_MEAN = (R0 + R1) / 2

        R_reluc[:, 0] = (np.log(R1 / R_MEAN)) / (self.mu0 * theta * self.la)
        R_reluc[:, 1] = 0.5 * theta / (self.mu0 * self.la * np.log(R1 / R0))
        R_reluc[:, 2] = (np.log(R_MEAN / R0)) / (self.mu0 * theta * self.la)
        R_reluc[:, 3] = 0.5 * theta / (self.mu0 * self.la * np.log(R1 / R0))
    else:
        print("Wrong mode, choice between ")

    return R_reluc
