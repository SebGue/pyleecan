# -*- coding: utf-8 -*-
import time
import numpy as np
from scipy.sparse import coo_matrix


def assembly_one_area(
    self,
    area_selected,
    reluc_list,
    list_elements_materials,
    Num_Unknowns,
    list_elem,
    permeability_element,
    boundary_condition_list,
):
    """Matrix massembler for one area in the motor (rotor or stator or airgap)

    Parameters
    ----------
    area_selected : integer
        area or materials selected to assembly.
    reluc_list : nd-array (size: m * 4 of floats)
        Values of each reluctance per cell
    list_geomerty: nd-array (size: m of integers)
        Material of each cell of the geometry
    Num_Unknowns : nd-array (size: n of integers)
        Unknowns in the linear system
    list_elem : nd-array (size: m of integers)
        Mesh cells in the motor geometry
    permeability_element : nd-array (size: m of floats)
        Permeability values of each cell
    boundary_condition_list : nd-array (size: n of integers)
        Data-strurture for boundary conditions evaluation

    Returns
    -------
    graph : Scipy CSR matrix
        sparse representation of the matrix.

    """
    # Assembly one area of the matrix
    N_total_element = self.Num_Unknowns.max() + 1

    mask_area = area_selected == list_elements_materials
    permeability_element = permeability_element[mask_area]
    reluc_list = reluc_list[mask_area]

    reluc_upper = reluc_list[:, 0] / permeability_element
    reluc_right = reluc_list[:, 1] / permeability_element
    reluc_down = reluc_list[:, 2] / permeability_element
    reluc_left = reluc_list[:, 3] / permeability_element

    # Select element with no boundary condition  and belonging "ara selected"
    index_unknown_1 = Num_Unknowns[list_elem[mask_area, 0]]
    mask_index_unknown_1 = index_unknown_1 != -1

    index_unknown_2 = Num_Unknowns[list_elem[mask_area, 1]]
    mask_index_unknown_2 = index_unknown_2 != -1

    index_unknown_3 = Num_Unknowns[list_elem[mask_area, 2]]
    mask_index_unknown_3 = index_unknown_3 != -1

    index_unknown_4 = Num_Unknowns[list_elem[mask_area, 3]]
    mask_index_unknown_4 = index_unknown_4 != -1

    if np.all(boundary_condition_list[list_elem[mask_area, :]] != 3):

        mask = mask_index_unknown_1 * mask_index_unknown_2
        I = index_unknown_1[mask]
        J = index_unknown_2[mask]
        data = -reluc_down[mask]
        graph = coo_matrix((data, (I, J)), shape=(N_total_element, N_total_element))

        mask = mask_index_unknown_2 * mask_index_unknown_3
        I = index_unknown_2[mask]
        J = index_unknown_3[mask]
        data = -reluc_right[mask]
        graph += coo_matrix((data, (I, J)), shape=(N_total_element, N_total_element))

        mask = mask_index_unknown_3 * mask_index_unknown_4
        I = index_unknown_3[mask]
        J = index_unknown_4[mask]
        data = -reluc_upper[mask]
        graph += coo_matrix((data, (I, J)), shape=(N_total_element, N_total_element))

        mask = mask_index_unknown_4 * mask_index_unknown_1
        I = index_unknown_4[mask]
        J = index_unknown_1[mask]
        data = -reluc_left[mask]
        graph += coo_matrix((data, (I, J)), shape=(N_total_element, N_total_element))
    else:

        mask_BC1 = boundary_condition_list[list_elem[mask_area, 0]] == 3
        mask_BC2 = boundary_condition_list[list_elem[mask_area, 1]] == 3
        mask_BC3 = boundary_condition_list[list_elem[mask_area, 2]] == 3
        mask_BC4 = boundary_condition_list[list_elem[mask_area, 3]] == 3

        plus_minus = np.array([-1, 1])

        mask = mask_index_unknown_1 * mask_index_unknown_2
        I = index_unknown_1[mask]
        J = index_unknown_2[mask]
        data = (
            reluc_down[mask]
            * plus_minus[(np.logical_xor(mask_BC1, mask_BC2)).astype(np.int32)][mask]
        )
        graph = coo_matrix((data, (I, J)), shape=(N_total_element, N_total_element))

        mask = mask_index_unknown_2 * mask_index_unknown_3
        I = index_unknown_2[mask]
        J = index_unknown_3[mask]
        data = (
            reluc_right[mask]
            * plus_minus[(np.logical_xor(mask_BC2, mask_BC3)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(N_total_element, N_total_element))

        mask = mask_index_unknown_3 * mask_index_unknown_4
        I = index_unknown_3[mask]
        J = index_unknown_4[mask]
        data = (
            reluc_upper[mask]
            * plus_minus[(np.logical_xor(mask_BC3, mask_BC4)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(N_total_element, N_total_element))

        mask = mask_index_unknown_4 * mask_index_unknown_1
        I = index_unknown_4[mask]
        J = index_unknown_1[mask]
        data = (
            reluc_left[mask]
            * plus_minus[(np.logical_xor(mask_BC4, mask_BC1)).astype(np.int32)][mask]
        )
        graph += coo_matrix((data, (I, J)), shape=(N_total_element, N_total_element))

    # add the symetric part
    graph += graph.T

    # Add diagonal
    index_unknown_1 = index_unknown_1[mask_index_unknown_1]
    index_unknown_2 = index_unknown_2[mask_index_unknown_2]
    index_unknown_3 = index_unknown_3[mask_index_unknown_3]
    index_unknown_4 = index_unknown_4[mask_index_unknown_4]

    I_Diag = np.unique(
        np.concatenate(
            [index_unknown_1, index_unknown_2, index_unknown_3, index_unknown_4]
        )
    )
    Diag = np.zeros(I_Diag.size)

    Diag[np.searchsorted(I_Diag, index_unknown_1)] = (
        reluc_down[mask_index_unknown_1] + reluc_left[mask_index_unknown_1]
    )
    Diag[np.searchsorted(I_Diag, index_unknown_2)] += (
        reluc_right[mask_index_unknown_2] + reluc_down[mask_index_unknown_2]
    )
    Diag[np.searchsorted(I_Diag, index_unknown_3)] += (
        reluc_upper[mask_index_unknown_3] + reluc_right[mask_index_unknown_3]
    )
    Diag[np.searchsorted(I_Diag, index_unknown_4)] += (
        reluc_upper[mask_index_unknown_4] + reluc_left[mask_index_unknown_4]
    )

    # Diag[np.in1d(I_Diag,index_unknown_1,assume_unique=True)] = reluc_down[mask_index_unknown_1]+reluc_left[mask_index_unknown_1]
    # Diag[np.in1d(I_Diag,index_unknown_2,assume_unique=True)] += reluc_right[mask_index_unknown_2]+reluc_down[mask_index_unknown_2]
    # Diag[np.in1d(I_Diag,index_unknown_3,assume_unique=True)] += reluc_upper[mask_index_unknown_3]+reluc_right[mask_index_unknown_3]
    # Diag[np.in1d(I_Diag,index_unknown_4,assume_unique=True)] += reluc_upper[mask_index_unknown_4]+reluc_left[mask_index_unknown_4]

    graph += coo_matrix(
        (Diag, (I_Diag, I_Diag)), shape=(N_total_element, N_total_element)
    )

    return graph
