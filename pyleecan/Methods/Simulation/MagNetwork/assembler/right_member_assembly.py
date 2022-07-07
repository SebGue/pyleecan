# -*- coding: utf-8 -*-
import time
import numpy as np
from scipy.sparse import coo_matrix


def right_member_assembly(
    self,
    cells_materials,
    Num_Unknowns,
    list_elem,
    list_coord,
    reluc_list,
    list_permeability,
    Br,
    mu0,
    la,
    mode,
    JA=None,
    JB=None,
    JC=None,
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

    mask_magnet = cells_materials == 4
    mu_pm = mu0 * list_permeability[4 - 1]
    nn = Num_Unknowns.max() + 1

    # np.savetxt("cell_materials.csv", cells_materials)
    # print("non-zero in mak magnet?", mask_magnet.sum())
    if mode == "cartesian":
        h_x = np.linalg.norm(
            list_coord[list_elem[:, 0]] - list_coord[list_elem[:, 1]], axis=1, ord=2
        )
        h_y = np.linalg.norm(
            list_coord[list_elem[:, 1]] - list_coord[list_elem[:, 2]], axis=1, ord=2
        )
        FMMPM = 0.5 * Br * la * h_y[mask_magnet] / mu_pm

    elif mode == "polar":
        R1 = list_coord[list_elem[:, 1], 0]
        R2 = list_coord[list_elem[:, 1], -1]
        sigma_R = np.abs(R2 - R1)
        sigma_theta = np.abs(
            list_coord[list_elem[:, 0], 1] - list_coord[list_elem[:, 0], 0]
        )

        FMMPM = 0.5 * Br * la * sigma_R[mask_magnet] / mu_pm

    ##Initailyze returned vector -> RHS
    E = np.zeros(nn, dtype=np.float64)

    ####Permanant Magnet assembling
    i1 = Num_Unknowns[list_elem[mask_magnet, 0]]
    i2 = Num_Unknowns[list_elem[mask_magnet, 1]]
    i3 = Num_Unknowns[list_elem[mask_magnet, 2]]
    i4 = Num_Unknowns[list_elem[mask_magnet, 3]]

    E[i1] += FMMPM * reluc_list[mask_magnet, 0]
    E[i2] -= FMMPM * reluc_list[mask_magnet, 1]
    E[i3] -= FMMPM * reluc_list[mask_magnet, 2]
    E[i4] += FMMPM * reluc_list[mask_magnet, 3]

    ####Winding assembling
    if JA is None and JB is None and JC is None:
        # Phase 1
        mask_winding = cells_materials == 1

        if mode == "cartesian":
            S = h_x[mask_winding] * h_y[mask_winding] / 4
        elif mode == "polar":
            S = (
                0.5
                * sigma_theta[mask_winding]
                * (R2[mask_winding] ** 2 - R1[mask_winding] ** 2)
                / 4
            )

        i1 = Num_Unknowns[list_elem[mask_winding, 0]]
        i2 = Num_Unknowns[list_elem[mask_winding, 1]]
        i3 = Num_Unknowns[list_elem[mask_winding, 2]]
        i4 = Num_Unknowns[list_elem[mask_winding, 3]]

        E[i1] += JA * S
        E[i2] += JA * S
        E[i3] += JA * S
        E[i4] += JA * S

        # Phase 2
        mask_winding = cells_materials == 2

        if mode == "cartesian":
            S = h_x[mask_winding] * h_y[mask_winding] / 4
        elif mode == "polar":
            S = (
                0.5
                * sigma_theta[mask_winding]
                * (R2[mask_winding] ** 2 - R1[mask_winding] ** 2)
                / 4
            )

        i1 = Num_Unknowns[list_elem[mask_winding, 0]]
        i2 = Num_Unknowns[list_elem[mask_winding, 1]]
        i3 = Num_Unknowns[list_elem[mask_winding, 2]]
        i4 = Num_Unknowns[list_elem[mask_winding, 3]]

        E[i1] += JB * S
        E[i2] += JB * S
        E[i3] += JB * S
        E[i4] += JB * S

        ###Phase 3
        mask_winding = cells_materials == 3
        if mode == "cartesian":
            S = h_x[mask_winding] * h_y[mask_winding] / 4
        elif mode == "polar":
            S = (
                0.5
                * sigma_theta[mask_winding]
                * (R2[mask_winding] ** 2 - R1[mask_winding] ** 2)
                / 4
            )

        i1 = Num_Unknowns[list_elem[mask_winding, 0]]
        i2 = Num_Unknowns[list_elem[mask_winding, 1]]
        i3 = Num_Unknowns[list_elem[mask_winding, 2]]
        i4 = Num_Unknowns[list_elem[mask_winding, 3]]

        E[i1] += JC * S
        E[i2] += JC * S
        E[i3] += JC * S
        E[i4] += JC * S

    return E