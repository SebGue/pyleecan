# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def plot_flux_line(self, Bx, By, cell_materials, m, n, x, y):
    """Plots the flux lines of the studied electric machine

    Parameters
    ----------
    self : Visualization
        A Visualization object
    Bx : float
        Magnetic flux density value in the x-direction (Default value = 0) [T]
    By : float
        Magnetic flux density value in the y-direction (Default value = 0) [T]
    cell_materials : int
        Cells used for the domain discretization (Default value = 1)
    n : int
        Number of discretization points along the x-axis (Default value = 1)
    m : int
        Number of discretization points along the y-axis (Default value = 1)
    x : float
        Coordinate of the point along the x-axis (Default value = 0) [m]
    y : float
        Coordinate of the point along the y-axis (Default value = 0) [m]

    Returns
    -------
    The flux lines mapping
        Flux Phi [wb]

    """

    "Mesh grid plotting"
    plt.pcolormesh(
        x,
        y,
        self.cell_materials.reshape((self.m, self.n)),
        shading="auto",
        cmap="Paired",
        edgecolors=None,
        facecolors="none",
        alpha=0.6,
    )

    "Flux contour plotting"
    ct = plt.contour(
        x,
        y,
        self.F.reshape((self.m + 1, self.n + 1)),
        cmap="jet",
        levels=18,
        alpha=1,
        linewidths=1,
    )

    "Figure title and legend"
    plt.title("Flux Lines Mapping")
    cbar = plt.colorbar(ct)
    cbar.ax.set_title("Flux $\phi (Weber)$")
    matplotlib.pyplot.xlabel("x-axis (m)")
    matplotlib.pyplot.ylabel("y-axis (m)")

    plt.show()

    return