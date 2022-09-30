# Method to compute the air gap magnetic flux density localy
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np


def comp_flux_airgap_local(
    self, r, theta, Phi, list_elem, list_coord, la, Rgap, type_coord_sys
):
    """
    Computes the flux density in the airgap

    Parameters
    ----------
    self : MagNetwork
        A MagNetwork object
    r: array
        Values of r coordinates
    theta: array
        Values of theta coordinates
    Phi: array
        Values of the flux Phi
    list_elem : ndarray
        List of the mesh cells inside the geometry
    list_coord : nd-array
        List of the geometry coordinate points
    la : integer
        Active length of the motor
    Rgap : integer
        Radius in the middle of the airgap
    type_coord_sys : integer (Default = 2)
        Type of the coordinate system : 1 for cartesian, 2 for Radial

    Returns
    -------
    Bx : nd-array
        Airgap flux density in the radial direction
    By: nd-array
        Airgap flux density in the theta direction
    """
    # Computing the radial and tangential flux density
    Bx, By = self.compute_B(Phi, list_elem, list_coord, la, type_coord_sys)

    # Reshaping B in a 2D array
    # Bx = Bx.reshape((r.size - 1, theta.size - 1))
    # By = By.reshape((r.size - 1, theta.size - 1))

    Bx = Bx.reshape((r.size, theta.size))
    By = By.reshape((r.size, theta.size))

    # Looking for the position in the center of the mecanical airgap
    position = Rgap

    # Getting the index of the line i from the position in the center of the mecanical airgap
    index_airgap = (int)((position - r.min()) * ((r.size - 1) / (r.max() - r.min())))+1

    return Bx[index_airgap, :], By[index_airgap, :]
