# -*- coding: utf-8 -*-

import numpy as np


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
    Num_Unknowns = -np.ones(BC_list.size, dtype=np.int32)
    mask = BC_list == 0
    mask[Periodic_point[:, 0]] = True

    Num_Unknowns[mask] = np.arange(mask.sum(), dtype=np.int32)
    Num_Unknowns[Periodic_point[:, 1]] = Num_Unknowns[Periodic_point[:, 0]]

    return Num_Unknowns
