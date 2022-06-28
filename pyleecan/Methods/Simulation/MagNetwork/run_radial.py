# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:46:17 2022

@author: LAP02
"""
# from solver_linear_model import linear_model
# import geometry_linear_motor
# from post_processing import compute_B_square
# from plot import view_contour_flux
import numpy as np
import meshio
import matplotlib.pyplot as plt
from pyleecan.Functions.load import load


def run_radial(self, axes_dict, Is_val=None):

    Machine = self.parent.machine
    la = 1  # Active length (m)
    Br = 1.2  # PM remanent induction (residual induction) (T)
    mu0 = np.pi * 4e-7  # Permeability of vacuum (H/m)

    pos = 2
    x_min = 0
    x_max = 0.06

    y_min = 0
    y_max = 0.051

    density = 1
    size_r = density * 51 + 1
    size_theta = density * 60 + 1

    # r = np.linspace(0.005, 0.05, size_r)
    # theta = np.linspace(0, np.pi / 2, size_theta)
    r = np.linspace(Machine.rotor.Rint, Machine.stator.Rext, size_r)
    # Add one extra point so that dual mesh has the correct dimension
    theta = np.linspace(
        axes_dict["angle"].initial,
        axes_dict["angle"].final,
        axes_dict["angle"].number + 1,
        endpoint=False,
    )

    r_dual = (r[1:] + r[:-1]) / 2
    theta_dual = (theta[1:] + theta[-1]) / 2

    BC = ["AP", "HD", "AP", "HD"]
    mode = "polar"

    # Compute current densities
    if Is_val is not None:
        surface_active = self.parent.machine.stator.slot.comp_surface_active()
        JA = Is_val[0, :] / surface_active
        JB = Is_val[1, :] / surface_active
        JC = Is_val[2, :] / surface_active
    else:
        JA = None
        JB = None
        JC = None

    (
        F,
        list_geometry,
        Num_unknowns,
        list_elem,
        permeability_cell,
        list_coord,
    ) = self.solver_linear_model(
        size_theta,
        size_r,
        theta,
        r,
        theta_dual,
        r_dual,
        pos,
        BC,
        self.geometry_linear_motor,
        mu0,
        la,
        Br,
        mode,
        JA=JA,
        JB=JB,
        JC=JC,
    )

    x = (list_coord[:, 1] * np.cos(list_coord[:, 0])).reshape(size_r, size_theta)
    y = (list_coord[:, 1] * np.sin(list_coord[:, 0])).reshape(size_r, size_theta)
    list_cart = np.stack((x, y), axis=-1)

    self.view_contour_flux(F, x, y, x.shape[1], x.shape[0], list_geometry)
    print(self.view_contour_flux(F, x, y, x.shape[1], x.shape[0], list_geometry))

    Bx, By = self.compute_B_square(F, list_elem, list_coord, la)

    B = np.stack((Bx, By), axis=-1)

    temp = np.zeros((list_elem.shape[0], 3))
    temp[:, 0] = Bx.flatten()
    temp[:, 1] = By.flatten()
    B = temp

    temp = np.zeros((list_coord.shape[0], 3))
    temp[:, 0:2] = list_coord

    list_coord = temp

    points = list_coord
    cells = [
        ("quad", list_elem),
    ]

    mesh = meshio.Mesh(
        points,
        cells,
        # Optionally provide extra data on points, cells, etc.
        point_data={"Flux": F},
        # Each item in cell data must match the cells array
        cell_data={"B": [B], "Materials": [list_geometry]},
    )
    mesh.write(
        "mymesh.xdmf",  # str, os.PathLike, or buffer/open file
        # file_format="vtk",  # optional if first argument is a path; inferred from extension
    )

    # Alternative with the same options
    # meshio.write_points_cells("mymesh.vtu", points, cells)

    print("mesh saved", list_coord.shape, list_elem.shape)

    # Compute 2D curve of the airgap flux density
    Bx_airgap, By_airgap = self.comp_flux_airgap_local(
        r, theta, F, list_elem, list_coord, la, Machine.stator.Rext, Machine.rotor.Rext
    )

    return Bx, By, Bx_airgap, By_airgap
