# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:02:37 2022

@author: LAP02
"""

import matplotlib.pyplot as plt
import numpy as np


def view_contour_flux2(self, F, x, y, size_x, size_y, cell_materials):
    # Plot flux
    plt.pcolormesh(
        x,
        y,
        cell_materials.reshape((size_y - 1, size_x - 1)),
        cmap="Paired",
        edgecolors=None,
        facecolors="none",
        alpha=0.6,
    )

    ct = plt.contour(
        F.reshape((size_y - 1, size_x - 1)),
        cmap="jet",
        levels=18,
        alpha=1,
        linewidths=1,
    )
    cbar = plt.colorbar(ct)
    cbar.ax.set_title("$Wb (Weber, Wb=T\cdot m^2)$")
    plt.show()
    return
