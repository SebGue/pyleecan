# -*- coding: utf-8 -*-
import numpy as np


def compute_B_radial(self, Phi, list_elem, list_coord, la):

    h_theta = np.abs(list_coord[list_elem[:, 0], 0] - list_coord[list_elem[:, 1], 0])

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
