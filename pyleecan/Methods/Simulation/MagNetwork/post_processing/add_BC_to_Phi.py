# -*- coding: utf-8 -*-
import numpy as np


def add_BC_to_Phi(self, Phi, Num_Unknowns, list_elem, boundary_condition_list):
    """
    Add boundary condition to the the flux.

    Parameters
    ----------
    Phi: nd-array, size: Num_Unknows.max()+1 (float)
        Flux
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    boundary_condition_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.

    Returns
    -------
    phi_total : size: n (float)
        The flux on each point of the mesh.

    """
    # post processing procedure to add boundary condition

    phi_total = np.zeros(Num_Unknowns.size)
    mask = Num_Unknowns != -1

    phi_total[mask] = Phi[Num_Unknowns[mask]]

    if 3 in boundary_condition_list:
        mask_AP = boundary_condition_list == 3
        phi_total[mask_AP] = -phi_total[mask_AP]

    return phi_total
