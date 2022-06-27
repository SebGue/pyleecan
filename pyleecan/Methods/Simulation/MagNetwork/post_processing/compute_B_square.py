# -*- coding: utf-8 -*-

import numpy as np


def compute_B_square(self, F, list_elem, list_coord, la):
    """
    Compute the magnetic field from the flux F and the mesh.

    Parameters
    ----------
    F : nd-array, size: n (float)
        Flux
    list_elem : nd-array, size: m (integers)
        tab of elements
    list_coord : nd-array, size: nx2 (flat)
        coordinate of each point.
    la : float
        depth of the model.

    Returns
    -------
    Bx : nd-array, size: m (flat)
        x component of B
    By : nd-array, size: m (float)
        y component of B

    """

    h_x = np.linalg.norm(
        list_coord[list_elem[:, 0]] - list_coord[list_elem[:, 1]],
        axis=1,
        ord=2,
    )
    h_y = np.linalg.norm(
        list_coord[list_elem[:, 1]] - list_coord[list_elem[:, 2]],
        axis=1,
        ord=2,
    )

    # Compute B
    Bx1 = (F[list_elem[:, 2]] - F[list_elem[:, 1]]) / (h_y * la)
    Bx2 = (F[list_elem[:, 3]] - F[list_elem[:, 0]]) / (h_y * la)

    By1 = (F[list_elem[:, 0]] - F[list_elem[:, 1]]) / (h_x * la)
    By2 = (F[list_elem[:, 3]] - F[list_elem[:, 2]]) / (h_x * la)

    Bx = (Bx1 + Bx2) / 2
    By = (By1 + By2) / 2

    return Bx, By
