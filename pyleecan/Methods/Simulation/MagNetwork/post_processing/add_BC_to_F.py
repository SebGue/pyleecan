# -*- coding: utf-8 -*-
import numpy as np


def add_BC_to_F(self, F, Num_Unknowns, list_elem, BC_list):
    """
    Add boundary condition to the the flux.

    Parameters
    ----------
    F : nd-array, size: Num_Unknows.max()+1 (float)
        Flux
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.

    Returns
    -------
    F_full : size: n (float)
        The flux on each point of the mesh.

    """
    # post processing procedure to add boundary condition

    F_full = np.zeros(Num_Unknowns.size)
    mask = Num_Unknowns != -1

    F_full[mask] = F[Num_Unknowns[mask]]

    if 3 in BC_list:
        mask_AP = BC_list == 3
        F_full[mask_AP] = -F_full[mask_AP]

    return F_full
