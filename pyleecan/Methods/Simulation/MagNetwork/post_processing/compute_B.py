# -*- coding: utf-8 -*-

import numpy as np


def compute_B(self, Phi, list_elem, list_coord, la, type_coord_sys):
    """
    Compute the magnetic field from the flux Phi and the mesh.

    Parameters
    ----------
    Phi : nd-array, size: n (float)
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
    if type_coord_sys == "cartesian":
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
        Bx_down = (Phi[list_elem[:, 2]] - Phi[list_elem[:, 1]]) / (h_y * la)
        Bx_upper = (Phi[list_elem[:, 3]] - Phi[list_elem[:, 0]]) / (h_y * la)

        By_left = (Phi[list_elem[:, 0]] - Phi[list_elem[:, 1]]) / (h_x * la)
        By_right = (Phi[list_elem[:, 3]] - Phi[list_elem[:, 2]]) / (h_x * la)

        Bx = (Bx_down + Bx_upper) / 2
        By = (By_left + By_right) / 2

        return Bx, By

    if type_coord_sys == "polar":
        h_theta = np.abs(
            list_coord[list_elem[:, 0], 0] - list_coord[list_elem[:, 1], 0]
        )
        r1 = list_coord[list_elem[:, 0], 1]
        r2 = list_coord[list_elem[:, 3], 1]
        h_r = np.abs(r2 - r1)

        # Compute B
        Btheta_down = (Phi[list_elem[:, 1]] - Phi[list_elem[:, 0]]) / (h_theta * r1 * la)
        Btheta_upper = (Phi[list_elem[:, 2]] - Phi[list_elem[:, 3]]) / (h_theta * r2 * la)

        Br_left = (Phi[list_elem[:, 0]] - Phi[list_elem[:, 3]]) / (h_r * la)
        Br_right = (Phi[list_elem[:, 1]] - Phi[list_elem[:, 2]]) / (h_r * la)

        Btheta = (Btheta_down + Btheta_upper) / 2
        Br = (Br_left + Br_right) / 2

        return Btheta, Br
