# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:02:37 2022

@author: LAP02
"""

import matplotlib.pyplot as plt
import numpy as np


def Plot_B(self, Bx, By, x, y, x_dual, y_dual, cell_materials, m, n):
    """

    Parameters
    ----------
    Bx : TYPE
        DESCRIPTION.
    By : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    x_dual : TYPE
        DESCRIPTION.
    y_dual : TYPE
        DESCRIPTION.
    cell_materials : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Getting the norm of the magnetix flux density B
    norm = np.sqrt(Bx ** 2 + By ** 2)

    k = 1
    skip = (slice(None, None, k), slice(None, None, k))

    # Defining the mesh grid (x, y)
    X, Y = np.meshgrid(x, y)

    # Plotting cell_materials
    plt.pcolormesh(
        X[skip],
        Y[skip],
        cell_materials.reshape((m, n))[skip],
        cmap="Paired",
        edgecolors=None,
        facecolors="none",
        alpha=0.6,
        linewidths=0.005,
    )

    # Defining the mesh grid (x_dual, y_dual)
    X, Y = np.meshgrid(x_dual, y_dual)

    # Plotting the magnetic flux density B
    plt.quiver(
        X[skip],
        Y[skip],
        Bx[skip] / norm[skip],
        By[skip] / norm[skip],
        norm[skip],
        cmap="cool",
        units="width",
        angles="xy",
        scale=30,
        width=0.01,
    )

    # Defining the colorbar
    cbar = plt.colorbar()

    # Defining the graph title
    cbar.ax.set_title("B (Tesla)")

    # Showing the final plot
    plt.show()
    return
