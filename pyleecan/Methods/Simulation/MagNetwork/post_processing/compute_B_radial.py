# -*- coding: utf-8 -*-
import numpy as np


def compute_B_radial(self, F, list_elem, list_coord, la):

    h_theta1 = np.abs(list_coord[list_elem[:, 0], 0] - list_coord[list_elem[:, 1], 0])
    h_theta2 = np.abs(list_coord[list_elem[:, 2], 0] - list_coord[list_elem[:, 3], 0])

    h_r = np.abs(list_coord[list_elem[:, 0], 1] - list_coord[list_elem[:, 3], 1])
    theta = list_coord[list_elem[:, 0], 1]

    # Compute B
    Btheta1 = (F[list_elem[:, 2]] - F[list_elem[:, 1]]) / (h_theta1 * theta * la)
    Btheta2 = (F[list_elem[:, 3]] - F[list_elem[:, 0]]) / (h_theta2 * theta * la)

    Br1 = (F[list_elem[:, 0]] - F[list_elem[:, 1]]) / (h_r * la)
    Br2 = (F[list_elem[:, 3]] - F[list_elem[:, 2]]) / (h_r * la)

    Btheta = (Btheta1 + Btheta2) / 2
    Br = (Br1 + Br2) / 2

    return Btheta, Br
