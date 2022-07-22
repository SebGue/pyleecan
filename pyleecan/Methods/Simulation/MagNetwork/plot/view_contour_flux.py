# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:02:37 2022

@author: LAP02
"""

import matplotlib.pyplot as plt
import numpy as np


def view_contour_flux(
    self, phi_matrix, theta, r, N_point_theta, N_point_r, list_element_materials
):

    try:
        # Plotting list_element_materials
        plt.pcolormesh(
            theta,
            r,
            list_element_materials.reshape((N_point_r - 1, N_point_theta - 1)),
            cmap="Paired",
            edgecolors=None,
            facecolors="none",
            alpha=0.6,
        )

        # plotting the flux phi_matrix
        ct = plt.contour(
            theta,
            r,
            phi_matrix.reshape((N_point_r, N_point_theta)),
            cmap="jet",
            levels=18,
            alpha=1,
            linewidths=1,
        )

        # Defining thecolorbar
        cbar = plt.colorbar(ct)

        # Defining the graph title
        cbar.ax.set_title("$Wb (Weber, Wb=T\cdot m^2)$")

        # Showing the final plot
        plt.show()

    except IndexError:
        print("Array is empty")
