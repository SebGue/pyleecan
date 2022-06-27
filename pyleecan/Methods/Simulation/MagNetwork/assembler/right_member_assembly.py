# -*- coding: utf-8 -*-
import time
import numpy as np
from scipy.sparse import coo_matrix


def right_member_assembly(
    self, cells_materials, Num_Unknowns, list_elem, list_coord, Br, mu0, mode
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
    mask_magnet = cells_materials == 4

    np.savetxt("cell_materials.csv", cells_materials)
    print("non-zero in mak magnet?", mask_magnet.sum())
    if mode == "cartesian":
        h_x = np.linalg.norm(
            list_coord[list_elem[:, 0]] - list_coord[list_elem[:, 1]], axis=1, ord=2
        )
        h_y = np.linalg.norm(
            list_coord[list_elem[:, 1]] - list_coord[list_elem[:, 2]], axis=1, ord=2
        )
        FMMPM = Br * h_x[0] * 0.5 / mu0

    elif mode == "polar":
        theta = np.abs(list_coord[list_elem[:, 0], 0] - list_coord[list_elem[:, 1], 0])
        R0 = list_coord[list_elem[:, 0], 1]
        R1 = list_coord[list_elem[:, 3], 1]

        h_x = 0.5 * (R0 + R1) * theta
        h_y = np.abs(R1 - R0)
        FMMPM = Br * np.mean(h_x[mask_magnet]) * 0.5 / mu0

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

    i1 = Num_Unknowns[list_elem[mask_magnet, 0]]
    i2 = Num_Unknowns[list_elem[mask_magnet, 1]]
    i3 = Num_Unknowns[list_elem[mask_magnet, 2]]
    i4 = Num_Unknowns[list_elem[mask_magnet, 3]]

    E[i1] += FMMPM
    E[i2] -= FMMPM
    E[i3] -= FMMPM
    E[i4] += FMMPM

    # mask_winding = cells_materials == 1
    # S = h_x[mask_winding] * h_y[mask_winding] / 4
    # i1 = Num_Unknowns[list_elem[mask_winding, 0]]
    # i2 = Num_Unknowns[list_elem[mask_winding, 1]]
    # i3 = Num_Unknowns[list_elem[mask_winding, 2]]
    # i4 = Num_Unknowns[list_elem[mask_winding, 3]]

    # E[i1] += JA * S
    # E[i2] += JA * S
    # E[i3] += JA * S
    # E[i4] += JA * S

    # mask_winding = cells_materials == 2
    # S = h_x[mask_winding] * h_y[mask_winding] / 4
    # i1 = Num_Unknowns[list_elem[mask_winding, 0]]
    # i2 = Num_Unknowns[list_elem[mask_winding, 1]]
    # i3 = Num_Unknowns[list_elem[mask_winding, 2]]
    # i4 = Num_Unknowns[list_elem[mask_winding, 3]]

    # E[i1] += JB * S
    # E[i2] += JB * S
    # E[i3] += JB * S
    # E[i4] += JB * S

    # mask_winding = cells_materials == 3
    # S = h_x[mask_winding] * h_y[mask_winding] / 4
    # i1 = Num_Unknowns[list_elem[mask_winding, 0]]
    # i2 = Num_Unknowns[list_elem[mask_winding, 1]]
    # i3 = Num_Unknowns[list_elem[mask_winding, 2]]
    # i4 = Num_Unknowns[list_elem[mask_winding, 3]]

    # E[i1] += JC * S
    # E[i2] += JC * S
    # E[i3] += JC * S
    # E[i4] += JC * S

    return E
