# -*- coding: utf-8 -*-

import numpy as np


def init_reluc(self, list_elem, list_coord, mu0, la, type_coord_sys):
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
    Reluc_list = np.zeros((list_elem.shape[0], 4))

    # Reluctances in the case of the cartesian coordiante system
    if type_coord_sys == 1:
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

        Reluc_list[:, 0] = 0.5 / (mu0 * la) * h_y / h_x
        Reluc_list[:, 1] = 0.5 / (mu0 * la) * h_x / h_y
        Reluc_list[:, 2] = 0.5 / (mu0 * la) * h_y / h_x
        Reluc_list[:, 3] = 0.5 / (mu0 * la) * h_x / h_y

    # Reluctances in the case of the polar coordiante system
    elif type_coord_sys == 2:
        theta = np.abs(list_coord[list_elem[:, 0], 0] - list_coord[list_elem[:, 1], 0])
        R0 = list_coord[list_elem[:, 0], 1]
        R1 = list_coord[list_elem[:, 3], 1]

        ln = np.log(R1 / R0)
        Reluc_list[:, 0] = 0.5 * ln / (mu0 * theta * la)
        Reluc_list[:, 1] = 0.5 * theta / (mu0 * la * ln)
        Reluc_list[:, 2] = Reluc_list[:, 0]
        Reluc_list[:, 3] = Reluc_list[:, 1]
    else:
        raise NameError("Wrong type_coord_sys, choice between 'cartesian' and 'polar'")

    return Reluc_list
