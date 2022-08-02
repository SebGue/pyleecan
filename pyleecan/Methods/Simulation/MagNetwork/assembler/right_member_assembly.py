# -*- coding: utf-8 -*-
import time
import numpy as np
from scipy.sparse import coo_matrix


def right_member_assembly(
    self,
    list_permeability_cell,  # cells materials
    Num_Unknowns,
    list_elem,
    list_coord,
    reluc_list,
    Br,
    mu0,
    la,
    type_coord_sys,
    JA=None,
    JB=None,
    JC=None,
):

    #######################################################################################
    # Number of PMs per period
    #######################################################################################
    # Get machine object
    Machine = self.parent.machine

    if Machine.comp_periodicity_spatial()[1] == True:
        periodicity = Machine.comp_periodicity_spatial()[0]
    else:
        periodicity = Machine.comp_periodicity_spatial()[0] / 2

    # Number of PMs per period
    nb_PM_per_period = round(Machine.rotor.get_pole_pair_number() / periodicity)

    # Relative permeabiltity of the PM
    mur_PM = Machine.rotor.magnet.mat_type.mag.mur_lin
    #######################################################################################

    mask_magnet = list_permeability_cell == mur_PM
    N_unknowns = Num_Unknowns.max() + 1

    #######################################################################################
    # Modeling of the PM force according to the PM direction (y): Phi_PM = FMMPM

    # Linear : Ref 1: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9473194
    if type_coord_sys == 1:
        h_x = np.linalg.norm(
            list_coord[list_elem[:, 0]] - list_coord[list_elem[:, 1]], axis=1, ord=2
        )
        h_y = np.linalg.norm(
            list_coord[list_elem[:, 1]] - list_coord[list_elem[:, 2]], axis=1, ord=2
        )
        FMMPM = 2 * Br * h_y[mask_magnet] * la / mur_PM

    # Radial: Ref 2: https://www.researchgate.net/publication/269405270_Modeling_of_a_Radial_Flux_PM_Rotating_Machine_using_a_New_Hybrid_Analytical_Model
    elif type_coord_sys == 2:
        R1 = list_coord[list_elem[:, 0], 1]
        R2 = list_coord[list_elem[:, -1], 1]
        sigma_R = np.abs(R2 - R1)
        sigma_theta = np.abs(
            list_coord[list_elem[:, 1], 0] - list_coord[list_elem[:, 0], 0]
        )

        FMMPM = 2 * Br * sigma_R[mask_magnet] * la / mur_PM

        # Initialize returned vector -> RHS
        RHS = np.zeros(N_unknowns, dtype=np.float64)

        # Permanant Magnet assembling
        index_unknowns_1 = Num_Unknowns[list_elem[mask_magnet, 0]]
        index_unknowns_2 = Num_Unknowns[list_elem[mask_magnet, 1]]
        index_unknowns_3 = Num_Unknowns[list_elem[mask_magnet, 2]]
        index_unknowns_4 = Num_Unknowns[list_elem[mask_magnet, 3]]

        RHS[index_unknowns_1] += FMMPM * reluc_list[mask_magnet, 3]
        RHS[index_unknowns_2] -= FMMPM * reluc_list[mask_magnet, 1]
        RHS[index_unknowns_3] -= FMMPM * reluc_list[mask_magnet, 1]
        RHS[index_unknowns_4] += FMMPM * reluc_list[mask_magnet, 3]

    ####Winding assembling

    if Machine.stator.winding.conductor.cond_mat.mag != None:
        mu_winding = Machine.stator.winding.conductor.cond_mat.mag.mur_lin
    else:
        mu_winding = 1

    if JA is None and JB is None and JC is None:
        # Phase 1
        mask_winding = list_permeability_cell == mu_winding

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
        mask_winding = list_permeability_cell == mu_winding + 1
        # here to be changed by the winding number bob1, bob2..

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
        mask_winding = list_permeability_cell == mu_winding + 2
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
