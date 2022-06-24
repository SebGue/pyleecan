# -*- coding: utf-8 -*-
import time
import numpy as np
from scipy.sparse import coo_matrix


def assembly(self, reluc_list, Num_Unknowns, list_elem, permeability_cell, BC_list):
    """
    Matrix assembler for all parts of the model

    Parameters
    ----------
    reluc_list : nd-array, size: mx4 (float)
        list contains the values of each reluctance per cell
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system
    list_elem : nd-array, size: m (integers)
        tab of elements
    permeability_cell : nd-array, size: m (float)
        the permeability values of each cell
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation

    Returns
    -------
    graph : Scipy CSR matrix
        sparse representation of the matrix

    """

    nn = Num_Unknowns.max() + 1

    R_N = reluc_list[:, 0] / permeability_cell
    R_E = reluc_list[:, 1] / permeability_cell
    R_S = reluc_list[:, 2] / permeability_cell
    R_W = reluc_list[:, 3] / permeability_cell

    # Select indices of each vertices of this elements (square)
    i1 = Num_Unknowns[list_elem[:, 0]]
    mask_i1 = i1 != -1
    i2 = Num_Unknowns[list_elem[:, 1]]
    mask_i2 = i2 != -1
    i3 = Num_Unknowns[list_elem[:, 2]]
    mask_i3 = i3 != -1
    i4 = Num_Unknowns[list_elem[:, 3]]
    mask_i4 = i4 != -1

    t1 = time.perf_counter()

    # Add tirangular part
    if np.all(BC_list != 3):

        mask = mask_i1 * mask_i2
        I = i1[mask]
        J = i2[mask]
        data = -R_S[mask]
        graph = coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i2 * mask_i3
        I = i2[mask]
        J = i3[mask]
        data = -R_E[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i3 * mask_i4
        I = i3[mask]
        J = i4[mask]
        data = -R_N[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i4 * mask_i1
        I = i4[mask]
        J = i1[mask]
        data = -R_W[mask]
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))
    else:
        mask_BC1 = BC_list[list_elem[:, 0]] == 3
        mask_BC2 = BC_list[list_elem[:, 1]] == 3
        mask_BC3 = BC_list[list_elem[:, 2]] == 3
        mask_BC4 = BC_list[list_elem[:, 3]] == 3

        plus_minus = np.array([-1, 1])

        mask = mask_i1 * mask_i2
        I = i1[mask]
        J = i2[mask]
        data = (
            R_S[mask]
            * plus_minus[(np.logical_xor(mask_BC1, mask_BC2)).astype(np.int32)][mask]
        )
        graph = coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i2 * mask_i3
        I = i2[mask]
        J = i3[mask]
        data = (
            R_E[mask]
            * plus_minus[(np.logical_xor(mask_BC2, mask_BC3)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i3 * mask_i4
        I = i3[mask]
        J = i4[mask]
        data = (
            R_N[mask]
            * plus_minus[(np.logical_xor(mask_BC3, mask_BC4)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

        mask = mask_i4 * mask_i1
        I = i4[mask]
        J = i1[mask]
        data = (
            R_W[mask]
            * plus_minus[(np.logical_xor(mask_BC4, mask_BC1)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(nn, nn))

    t2 = time.perf_counter()
    print("time assembly matrix without diag", t2 - t1)
    # add the symetric part
    graph += graph.T

    # Add diagonal
    Diag = np.zeros(nn, dtype=np.float64)

    Diag[i1[mask_i1]] += R_S[mask_i1] + R_W[mask_i1]
    Diag[i2[mask_i2]] += R_E[mask_i2] + R_S[mask_i2]
    Diag[i3[mask_i3]] += R_N[mask_i3] + R_E[mask_i3]
    Diag[i4[mask_i4]] += R_N[mask_i4] + R_W[mask_i4]

    I = np.arange(nn)
    graph += coo_matrix((Diag, (I, I)), shape=(nn, nn))
    # graph.setdiag(Diag)
    t3 = time.perf_counter()
    print("time assembly diag", t3 - t2)

    return graph
