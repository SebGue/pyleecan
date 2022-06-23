# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:11:39 2022

@author: LAP02
"""

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
        self.list_coord[self.list_elem[:, 0]] - self.list_coord[self.list_elem[:, 1]],
        axis=1,
        ord=2,
    )
    h_y = np.linalg.norm(
        self.list_coord[self.list_elem[:, 1]] - self.list_coord[self.list_elem[:, 2]],
        axis=1,
        ord=2,
    )

    # Compute B
    Bx1 = (self.F[self.list_elem[:, 2]] - self.F[self.list_elem[:, 1]]) / (
        h_y * self.la
    )
    Bx2 = (self.F[self.list_elem[:, 3]] - self.F[self.list_elem[:, 0]]) / (
        h_y * self.la
    )

    By1 = (self.F[self.list_elem[:, 0]] - self.F[self.list_elem[:, 1]]) / (
        h_x * self.la
    )
    By2 = (self.F[self.list_elem[:, 3]] - self.F[self.list_elem[:, 2]]) / (
        h_x * self.la
    )

    Bx = (Bx1 + Bx2) / 2
    By = (By1 + By2) / 2

    return Bx, By


def compute_B_radial(self, F, list_elem, list_coord, la):

    h_theta1 = np.abs(
        self.list_coord[self.list_elem[:, 0], 0]
        - self.list_coord[self.list_elem[:, 1], 0]
    )
    h_theta2 = np.abs(
        self.list_coord[self.list_elem[:, 2], 0]
        - self.list_coord[self.list_elem[:, 3], 0]
    )

    h_r = np.abs(
        self.list_coord[self.list_elem[:, 0], 1]
        - self.list_coord[self.list_elem[:, 3], 1]
    )
    theta = self.list_coord[self.list_elem[:, 0], 1]

    # Compute B
    Btheta1 = (self.F[self.list_elem[:, 2]] - self.F[self.list_elem[:, 1]]) / (
        h_theta1 * theta * self.la
    )
    Btheta2 = (self.F[self.list_elem[:, 3]] - self.F[self.list_elem[:, 0]]) / (
        h_theta2 * theta * self.la
    )

    Br1 = (self.F[self.list_elem[:, 0]] - self.F[self.list_elem[:, 1]]) / (
        h_r * self.la
    )
    Br2 = (self.F[self.list_elem[:, 3]] - self.F[self.list_elem[:, 2]]) / (
        h_r * self.la
    )

    Btheta = (Btheta1 + Btheta2) / 2
    Br = (Br1 + Br2) / 2

    return Btheta, Br


def add_BC_to_F(self, F, Num_Unknowns, list_elem, BC_list):
    """
    Add boundary condition to the the flux.

    Parameters
    ----------
    F : nd-array, size: Num_Unknows.max()+1 (float)
        Flux
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.

    Returns
    -------
    F_full : size: n (float)
        The flux on each point of the mesh.

    """
    # post processing procedure to add boundary condition

    F_full = np.zeros(self.Num_Unknowns.size)
    mask = self.Num_Unknowns != -1

    F_full[mask] = self.F[self.Num_Unknowns[mask]]

    if 3 in self.BC_list:
        mask_AP = self.BC_list == 3
        F_full[mask_AP] = -F_full[mask_AP]

    return F_full
