# -*- coding: utf-8 -*-

import numpy as np


def init_permeabilty_cell(
    self, N_point_theta, N_point_r, permeability_materials, mu0, list_element_materials
):
    """
    Initaialyze the permeabilty in each cell.
    Parameters
    ----------
    N_point_theta : integer
        Number of horizontal points.
    N_point_r : integer
        Number of vertical point.
    permeability_materials : nd-array, size: (foat)
        permabilty of each materials.
    list_element_materials : nd-array, size: m (int)
        materials in each cell.
    Returns
    -------
    list_element_materials : nd-array, size: n (int)
        Permeability in each cell.
    """
    list_elem_permeability = np.zeros(
        (N_point_r - 1) * (N_point_theta - 1), dtype=np.float64
    )
    list_elem_permeability[...] = permeability_materials[list_element_materials - 1]

    return list_elem_permeability
