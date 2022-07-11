# -*- coding: utf-8 -*-

import numpy as np


def init_permeabilty_cell(
    self, size_x, size_y, permeability_materials, mu0, cells_materials
):
    """
    Initaialyze the permeabilty in each cell.

    Parameters
    ----------
    size_x : integer
        Number of horizontal points.
    size_y : integer
        Number of vertical point.
    permeability_materials : nd-array, size: (foat)
        permabilty of each materials.
    cells_materials : nd-array, size: m (int)
        materials in each cell.

    Returns
    -------
    permeability_cell : nd-array, size: n (int)
        Permeability in each cell.

    """
    permeability_cell = np.zeros((size_y - 1) * (size_x - 1), dtype=np.float64)
    # permeability_cell[...] = permeability_materials[cells_materials - 1] * mu0
    for i in range(len(permeability_cell)):
        permeability_cell[i] = permeability_materials[cells_materials[i] - 1] * mu0

    return permeability_cell
