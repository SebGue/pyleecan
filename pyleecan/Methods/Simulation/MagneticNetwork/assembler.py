# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:20:12 2022

@author: LAP02
"""

import time
import numpy as np
from scipy.sparse import coo_matrix


def assembly(reluc_list, Num_Unknowns, list_elem, permeability_cell, BC_list):
    """
    Matrix massembler for all part of the model

    Parameters
    ----------
    reluc_list : nd-array, size: mx4 (float)
        list contains the values of each reluctance by cell.
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    permeability_cell : nd-array, size: m (float)
        the permeability values of each cell
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.


    Returns
    -------
    graph : Scipy CSR matrix
        sparse representation of the matrix.

    """

    nn = Num_Unknowns.max() + 1

    R_N = reluc_list[:, 0] / permeability_cell
    R_E = reluc_list[:, 1] / permeability_cell
    R_S = reluc_list[:, 2] / permeability_cell
    R_W = reluc_list[:, 3] / permeability_cell

    # Select indices of each vertices of this elements (square)
    i1 = Num_Unknowns[list_elem[:, 0]]
    mask_i1 = i1 != -1
    i2 = Num_Unknowns[list_elem[:, 1]]
    mask_i2 = i2 != -1
    i3 = Num_Unknowns[list_elem[:, 2]]
    mask_i3 = i3 != -1
    i4 = Num_Unknowns[list_elem[:, 3]]
    mask_i4 = i4 != -1

    t1 = time.perf_counter()

    # Add tirangular part
    if np.all(BC_list != 3):

        mask = mask_i1 * mask_i2
        I = i1[mask]
        J = i2[mask]
        data = -R_S[mask]
        graph = coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i2 * mask_i3
        I = i2[mask]
        J = i3[mask]
        data = -R_E[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i3 * mask_i4
        I = i3[mask]
        J = i4[mask]
        data = -R_N[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i4 * mask_i1
        I = i4[mask]
        J = i1[mask]
        data = -R_W[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))
    else:
        mask_BC1 = BC_list[list_elem[:, 0]] == 3
        mask_BC2 = BC_list[list_elem[:, 1]] == 3
        mask_BC3 = BC_list[list_elem[:, 2]] == 3
        mask_BC4 = BC_list[list_elem[:, 3]] == 3

        plus_minus = np.array([-1, 1])

        mask = mask_i1 * mask_i2
        I = i1[mask]
        J = i2[mask]
        data = (
            R_S[mask]
            * plus_minus[(np.logical_xor(mask_BC1, mask_BC2)).astype(np.int32)][mask]
        )
        graph = coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i2 * mask_i3
        I = i2[mask]
        J = i3[mask]
        data = (
            R_E[mask]
            * plus_minus[(np.logical_xor(mask_BC2, mask_BC3)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i3 * mask_i4
        I = i3[mask]
        J = i4[mask]
        data = (
            R_N[mask]
            * plus_minus[(np.logical_xor(mask_BC3, mask_BC4)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i4 * mask_i1
        I = i4[mask]
        J = i1[mask]
        data = (
            R_W[mask]
            * plus_minus[(np.logical_xor(mask_BC4, mask_BC1)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

    t2 = time.perf_counter()
    print("time assembly matrix without diag", t2 - t1)
    # add the symetric part
    graph += graph.T

    # Add diagonal
    Diag = np.zeros(nn, dtype=np.float64)

    Diag[i1[mask_i1]] += R_S[mask_i1] + R_W[mask_i1]
    Diag[i2[mask_i2]] += R_E[mask_i2] + R_S[mask_i2]
    Diag[i3[mask_i3]] += R_N[mask_i3] + R_E[mask_i3]
    Diag[i4[mask_i4]] += R_N[mask_i4] + R_W[mask_i4]

    I = np.arange(nn)
    graph += coo_matrix((Diag, (I, I)), shape=(nn, nn))
    # graph.setdiag(Diag)
    t3 = time.perf_counter()
    print("time assembly diag", t3 - t2)

    return graph


def assembly_one_area(
    area_selected,
    reluc_list,
    list_geometry,
    Num_Unknowns,
    list_elem,
    permeability_cell,
    BC_list,
):
    """
    Matrix massembler for one area (materials) selected.

    Parameters
    ----------
    area_selected : integer
        area or materials selected to assembly.
    reluc_list : nd-array, size: mx4 (float)
        list contains the values of each reluctance by cell.
    list_geomerty: nd-array, size: m (integers)
        contains the material of each cells (elements)
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    permeability_cell : nd-array, size: m (float)
        the permeability values of each cell
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.


    Returns
    -------
    graph : Scipy CSR matrix
        sparse representation of the matrix.

    """
    # Assembly one area of the matrix
    nn = Num_Unknowns.max() + 1

    mask_area = area_selected == list_geometry
    permeability_cell = permeability_cell[mask_area]
    reluc_list = reluc_list[mask_area]

    R_N = reluc_list[:, 0] / permeability_cell
    R_E = reluc_list[:, 1] / permeability_cell
    R_S = reluc_list[:, 2] / permeability_cell
    R_W = reluc_list[:, 3] / permeability_cell

    # Select element with no boundary condition  and belonging "ara selected"
    i1 = Num_Unknowns[list_elem[mask_area, 0]]
    mask_i1 = i1 != -1

    i2 = Num_Unknowns[list_elem[mask_area, 1]]
    mask_i2 = i2 != -1

    i3 = Num_Unknowns[list_elem[mask_area, 2]]
    mask_i3 = i3 != -1

    i4 = Num_Unknowns[list_elem[mask_area, 3]]
    mask_i4 = i4 != -1

    if np.all(BC_list[list_elem[mask_area, :]] != 3):

        mask = mask_i1 * mask_i2
        I = i1[mask]
        J = i2[mask]
        data = -R_S[mask]
        graph = coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i2 * mask_i3
        I = i2[mask]
        J = i3[mask]
        data = -R_E[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i3 * mask_i4
        I = i3[mask]
        J = i4[mask]
        data = -R_N[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i4 * mask_i1
        I = i4[mask]
        J = i1[mask]
        data = -R_W[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))
    else:

        mask_BC1 = BC_list[list_elem[mask_area, 0]] == 3
        mask_BC2 = BC_list[list_elem[mask_area, 1]] == 3
        mask_BC3 = BC_list[list_elem[mask_area, 2]] == 3
        mask_BC4 = BC_list[list_elem[mask_area, 3]] == 3

        plus_minus = np.array([-1, 1])

        mask = mask_i1 * mask_i2
        I = i1[mask]
        J = i2[mask]
        data = (
            R_S[mask]
            * plus_minus[(np.logical_xor(mask_BC1, mask_BC2)).astype(np.int32)][mask]
        )
        graph = coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i2 * mask_i3
        I = i2[mask]
        J = i3[mask]
        data = (
            R_E[mask]
            * plus_minus[(np.logical_xor(mask_BC2, mask_BC3)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i3 * mask_i4
        I = i3[mask]
        J = i4[mask]
        data = (
            R_N[mask]
            * plus_minus[(np.logical_xor(mask_BC3, mask_BC4)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i4 * mask_i1
        I = i4[mask]
        J = i1[mask]
        data = (
            R_W[mask]
            * plus_minus[(np.logical_xor(mask_BC4, mask_BC1)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

    # add the symetric part
    graph += graph.T

    # Add diagonal
    i1 = i1[mask_i1]
    i2 = i2[mask_i2]
    i3 = i3[mask_i3]
    i4 = i4[mask_i4]

    I_Diag = np.unique(np.concatenate([i1, i2, i3, i4]))
    Diag = np.zeros(I_Diag.size)

    Diag[np.searchsorted(I_Diag, i1)] = R_S[mask_i1] + R_W[mask_i1]
    Diag[np.searchsorted(I_Diag, i2)] += R_E[mask_i2] + R_S[mask_i2]
    Diag[np.searchsorted(I_Diag, i3)] += R_N[mask_i3] + R_E[mask_i3]
    Diag[np.searchsorted(I_Diag, i4)] += R_N[mask_i4] + R_W[mask_i4]

    # Diag[np.in1d(I_Diag,i1,assume_unique=True)] = R_S[mask_i1]+R_W[mask_i1]
    # Diag[np.in1d(I_Diag,i2,assume_unique=True)] += R_E[mask_i2]+R_S[mask_i2]
    # Diag[np.in1d(I_Diag,i3,assume_unique=True)] += R_N[mask_i3]+R_E[mask_i3]
    # Diag[np.in1d(I_Diag,i4,assume_unique=True)] += R_N[mask_i4]+R_W[mask_i4]

    graph += coo_matrix((Diag, (I_Diag, I_Diag)), shape=(nn, nn))

    return graph


def right_member_assembly(
    cells_materials, Num_Unknowns, list_elem, list_coord, Br, mu0, mode
):
    """
    Vector assembler, RHS.

    Parameters
    ----------
    cells_materials: nd-array, size: m (integers)
        contains the material of each cells (elements)
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    list_coord : nd-array, size: nx2 (flat)
        coordinate of each point.
    Br : float
        Caracterization of the permanent magnet.
    mu0 : float
        Permaebility in vaccum

    Returns
    -------
    E : nd-array, size: Num.unknowns.max()+1
        RHS

    """
    nn = Num_Unknowns.max() + 1
    if mode == "cartesian":
        h_x = np.linalg.norm(
            list_coord[list_elem[:, 0]] - list_coord[list_elem[:, 1]], axis=1, ord=2
        )
        h_y = np.linalg.norm(
            list_coord[list_elem[:, 1]] - list_coord[list_elem[:, 2]], axis=1, ord=2
        )

    elif mode == "polar":
        theta = np.abs(list_coord[list_elem[:, 0], 0] - list_coord[list_elem[:, 1], 0])
        R0 = list_coord[list_elem[:, 0], 1]
        R1 = list_coord[list_elem[:, 3], 1]

        h_x = 0.5 * (R0 + R1) * theta
        h_y = np.abs(R1 - R0)
        print("ok", h_x, h_y)

    # FMMPM=0
    # print("FMMPM:", FMMPM)
    J = 5e6
    # J=0
    wt = 0
    E = np.zeros(nn, dtype=np.float64)
    JA = -J * np.cos(wt - 2 * np.pi / 3)
    JC = -J * np.cos(wt + 2 * np.pi / 3)
    JB = J * np.cos(wt)
    JA, JB, JC = 0, 0, 0

    mask_magnet = cells_materials == 4

    FMMPM = Br * h_x[0] * 0.5 / mu0
    # FMMPM = 0
    i1 = Num_Unknowns[list_elem[mask_magnet, 0]]
    i2 = Num_Unknowns[list_elem[mask_magnet, 1]]
    i3 = Num_Unknowns[list_elem[mask_magnet, 2]]
    i4 = Num_Unknowns[list_elem[mask_magnet, 3]]

    E[i1] -= FMMPM
    E[i2] += FMMPM
    E[i3] += FMMPM
    E[i4] -= FMMPM

    mask_winding = cells_materials == 1
    S = h_x[mask_winding] * h_y[mask_winding] / 4
    i1 = Num_Unknowns[list_elem[mask_winding, 0]]
    i2 = Num_Unknowns[list_elem[mask_winding, 1]]
    i3 = Num_Unknowns[list_elem[mask_winding, 2]]
    i4 = Num_Unknowns[list_elem[mask_winding, 3]]

    E[i1] += JA * S
    E[i2] += JA * S
    E[i3] += JA * S
    E[i4] += JA * S

    mask_winding = cells_materials == 2
    S = h_x[mask_winding] * h_y[mask_winding] / 4
    i1 = Num_Unknowns[list_elem[mask_winding, 0]]
    i2 = Num_Unknowns[list_elem[mask_winding, 1]]
    i3 = Num_Unknowns[list_elem[mask_winding, 2]]
    i4 = Num_Unknowns[list_elem[mask_winding, 3]]

    E[i1] += JB * S
    E[i2] += JB * S
    E[i3] += JB * S
    E[i4] += JB * S

    mask_winding = cells_materials == 3
    S = h_x[mask_winding] * h_y[mask_winding] / 4
    i1 = Num_Unknowns[list_elem[mask_winding, 0]]
    i2 = Num_Unknowns[list_elem[mask_winding, 1]]
    i3 = Num_Unknowns[list_elem[mask_winding, 2]]
    i4 = Num_Unknowns[list_elem[mask_winding, 3]]

    E[i1] += JC * S
    E[i2] += JC * S
    E[i3] += JC * S
    E[i4] += JC * S

    return E
