# -*- coding: utf-8 -*-

import numpy as np


def init_reluc(self, list_elem, list_coord, mu0, la, mode):
    """
    Parameters
    ----------
    list_elem : nd-array of integers ( m x 4 )
        list of the element
    list_coord : nd-array of coordinate ( n x 2 )
        coordinate of each vertices

    Returns
    -------
    Reluc_list : Reluc_list (m x 4 )
        Reluctance of each elements
    """
    R_reluc = np.zeros((list_elem.shape[0], 4))

    # Reluctances in the case of the cartesian coordiante system
    if mode == "cartesian":
        # length and width
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

        R_reluc[:, 0] = 0.5 / (mu0 * la) * h_y / h_x
        R_reluc[:, 1] = 0.5 / (mu0 * la) * h_x / h_y
        R_reluc[:, 2] = 0.5 / (mu0 * la) * h_y / h_x
        R_reluc[:, 3] = 0.5 / (mu0 * la) * h_x / h_y

    # Reluctances in the case of the polar coordiante system
    elif mode == "polar":
        theta = np.abs(list_coord[list_elem[:, 0], 0] - list_coord[list_elem[:, 1], 0])
        R0 = list_coord[list_elem[:, 0], 1]
        R1 = list_coord[list_elem[:, 3], 1]


        ln=np.log(R1 / R0)
        R_reluc[:, 0] = 0.5* ln / (mu0 * theta * la)
        R_reluc[:, 1] = 0.5 * theta / (mu0 * la * ln)
        R_reluc[:, 2] = R_reluc[:, 0]
        R_reluc[:, 3] = R_reluc[:, 1]
    else:
        raise NameError("Wrong mode, choice between 'cartesian' and 'polar'")

    return R_reluc
