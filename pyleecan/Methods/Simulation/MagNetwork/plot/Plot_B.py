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
    # Plot Magnetic field
    norm = np.sqrt(Bx ** 2 + By ** 2)

    k = 1
    skip = (slice(None, None, k), slice(None, None, k))
    X, Y = np.meshgrid(x, y)
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
    X, Y = np.meshgrid(x_dual, y_dual)
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
    cbar = plt.colorbar()
    cbar.ax.set_title("B (Tesla)")
    plt.show()
    return
