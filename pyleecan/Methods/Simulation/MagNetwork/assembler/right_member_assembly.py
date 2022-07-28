# -*- coding: utf-8 -*-
import time
import numpy as np
from scipy.sparse import coo_matrix


def right_member_assembly(
    self,
    list_elements_materials,
    Num_Unknowns,
    list_elem,
    list_coord,
    reluc_list,
    list_permeability,
    Br,
    mu0,
    la,
    type_coord_sys,
    JA=None,
    JB=None,
    JC=None,
):
    """Vector assembler (RHS)

    Parameters
    ----------
    list_elements_materials: nd-array (size: m of integers)
        Material of each cell of the motor geometry
    Num_Unknowns : nd-array (size: n of integers)
        Unknowns in the linear system
    list_elem : nd-array (size: m of integers)
        Elements in the motor geometry
    list_coord : nd-array (size: n * 2 (flat))
        Coordinates of each point
    Br : float
        Remeanance flux density of the permanent magnet
    mu0 : float
        Permeability of the vaccum

    Returns
    -------
    RHS : nd-array, size: Num.unknowns.max()+1
        Right member vector

    """

    position_permeability_PM = len(list_permeability) - 3
    mask_magnet = list_elements_materials == position_permeability_PM
    mur_PM = list_permeability[position_permeability_PM]  # 4-1???????
    N_unknowns = Num_Unknowns.max() + 1

    # np.savetxt("cell_materials.csv", list_elements_materials)
    # print("non-zero in mak magnet?", mask_magnet.sum())
    if type_coord_sys == 1:
        h_x = np.linalg.norm(
            list_coord[list_elem[:, 0]] - list_coord[list_elem[:, 1]], axis=1, ord=2
        )
        h_y = np.linalg.norm(
            list_coord[list_elem[:, 1]] - list_coord[list_elem[:, 2]], axis=1, ord=2
        )
        FMMPM = 0.5 * Br * la * h_y[mask_magnet] / mur_PM

    elif type_coord_sys == 2:
        R1 = list_coord[list_elem[:, 0], 1]
        R2 = list_coord[list_elem[:, -1], 1]
        sigma_R = np.abs(R2 - R1)
        sigma_theta = np.abs(
            list_coord[list_elem[:, 1], 0] - list_coord[list_elem[:, 0], 0]
        )

        FMMPM = 1.05 * Br * la * sigma_R[mask_magnet] * 2

        ##Initailyze returned vector -> RHS
        RHS = np.zeros(N_unknowns, dtype=np.float64)

        ####Permanant Magnet assembling
        index_unknowns_1 = Num_Unknowns[list_elem[mask_magnet, 0]]
        index_unknowns_2 = Num_Unknowns[list_elem[mask_magnet, 1]]
        index_unknowns_3 = Num_Unknowns[list_elem[mask_magnet, 2]]
        index_unknowns_4 = Num_Unknowns[list_elem[mask_magnet, 3]]

        RHS[index_unknowns_1] += FMMPM * reluc_list[mask_magnet, 3]
        RHS[index_unknowns_2] -= FMMPM * reluc_list[mask_magnet, 1]
        RHS[index_unknowns_3] -= FMMPM * reluc_list[mask_magnet, 1]
        RHS[index_unknowns_4] += FMMPM * reluc_list[mask_magnet, 3]

    ####Winding assembling
    if JA is None and JB is None and JC is None:
        # Phase 1
        mask_winding = list_elements_materials == 1

        if type_coord_sys == 1:
            S = h_x[mask_winding] * h_y[mask_winding] / 4
        elif type_coord_sys == 2:
            S = (
                0.5
                * sigma_theta[mask_winding]
                * (R2[mask_winding] ** 2 - R1[mask_winding] ** 2)
                / 4
            )

        index_unknowns_1 = Num_Unknowns[list_elem[mask_winding, 0]]
        index_unknowns_2 = Num_Unknowns[list_elem[mask_winding, 1]]
        index_unknowns_3 = Num_Unknowns[list_elem[mask_winding, 2]]
        index_unknowns_4 = Num_Unknowns[list_elem[mask_winding, 3]]

        RHS[index_unknowns_1] += JA * S
        RHS[index_unknowns_2] += JA * S
        RHS[index_unknowns_3] += JA * S
        RHS[index_unknowns_4] += JA * S

        # Phase 2
        mask_winding = list_elements_materials == 2

        if type_coord_sys == 1:
            S = h_x[mask_winding] * h_y[mask_winding] / 4
        elif type_coord_sys == 2:
            S = (
                0.5
                * sigma_theta[mask_winding]
                * (R2[mask_winding] ** 2 - R1[mask_winding] ** 2)
                / 4
            )

        index_unknowns_1 = Num_Unknowns[list_elem[mask_winding, 0]]
        index_unknowns_2 = Num_Unknowns[list_elem[mask_winding, 1]]
        index_unknowns_3 = Num_Unknowns[list_elem[mask_winding, 2]]
        index_unknowns_4 = Num_Unknowns[list_elem[mask_winding, 3]]

        RHS[index_unknowns_1] += JB * S
        RHS[index_unknowns_2] += JB * S
        RHS[index_unknowns_3] += JB * S
        RHS[index_unknowns_4] += JB * S

        ###Phase 3
        mask_winding = list_elements_materials == 3
        if type_coord_sys == 1:
            S = h_x[mask_winding] * h_y[mask_winding] / 4
        elif type_coord_sys == 2:
            S = (
                0.5
                * sigma_theta[mask_winding]
                * (R2[mask_winding] ** 2 - R1[mask_winding] ** 2)
                / 4
            )

        index_unknowns_1 = Num_Unknowns[list_elem[mask_winding, 0]]
        index_unknowns_2 = Num_Unknowns[list_elem[mask_winding, 1]]
        index_unknowns_3 = Num_Unknowns[list_elem[mask_winding, 2]]
        index_unknowns_4 = Num_Unknowns[list_elem[mask_winding, 3]]

        RHS[index_unknowns_1] += JC * S
        RHS[index_unknowns_2] += JC * S
        RHS[index_unknowns_3] += JC * S
        RHS[index_unknowns_4] += JC * S
    # np.savetxt("Everif.csv",RHS.reshape((50,60)),fmt="%5.2f")

    return RHS
