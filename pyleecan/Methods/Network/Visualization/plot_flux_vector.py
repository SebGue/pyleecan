# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def plot_flux_vector(self, F, m, n, cell_materials, x, y, x_dual, y_dual):
    """Plots the flux vectors of the studied electric machine

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
    x_dual : float
        Coordinate of center cell along the x-axis (Default value = 0) [m]
    y_dual : float
        Coordinate of center cell along the y-axis (Default value = 0) [m]

    Returns
    -------
    The flux density vectors mapping
        Flux Density [T]

    """

    """ Norm of the magnetic flux density B """
    norm = np.sqrt(self.Bx ** 2 + self.By ** 2)

    """ Slicing of the 2D plot for visible vectors plotting """
    k = 4
    skip = (slice(None, None, k), slice(None, None, k))

    """ Mesh grid plotting """
    X, Y = np.meshgrid(x, y)
    plt.pcolormesh(
        X[skip],
        Y[skip],
        self.cell_materials.reshape((self.m, self.n))[skip],
        shading="auto",
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
        self.Bx[skip] / norm[skip],
        self.By[skip] / norm[skip],
        norm[skip],
        cmap="cool",
        units="width",
        angles="xy",
        scale=30,
        width=0.006,
    )

    "Figure title and legend"
    plt.title("Flux Density Mapping")
    cbar = plt.colorbar()
    cbar.ax.set_title("Flux Density B (Tesla)")
    matplotlib.pyplot.xlabel("x-axis (m)")
    matplotlib.pyplot.ylabel("y-axis (m)")

    plt.show()
    return
