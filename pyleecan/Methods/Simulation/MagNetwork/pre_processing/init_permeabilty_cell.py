# -*- coding: utf-8 -*-

import numpy as np


def init_permeabilty_cell(
    self, N_point_theta, N_point_r, permeability_materials, mu0, list_element_materials
):
    """
    Initialize the permeabilties in each cell.

    Parameters
    ----------
    N_point_theta : integer
        Discretization point in the theta direction
    N_point_r : integer
        Discretization points in the r direction
    permeability_materials : nd-array of floats
        Permabilty of each material contituting the motor
    list_element_materials : nd-array (size: n of integers)
        Permeability in each cell inside the motor geometry

    Returns
    -------
    list_elem_permeability : nd-array (size: n of integers)
        Permeability in each cell inside the motor geometry
    """
    list_elem_permeability = np.zeros(
        (N_point_r - 1) * (N_point_theta - 1), dtype=np.float64
    )

    list_elem_permeability[...] = permeability_materials[list_element_materials - 1]

    return list_elem_permeability
