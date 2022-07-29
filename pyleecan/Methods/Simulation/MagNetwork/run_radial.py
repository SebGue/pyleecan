# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:46:17 2022
@author: LAP02
"""
from re import S
import numpy as np
import meshio
import matplotlib.pyplot as plt
from pyleecan.Functions.load import load

from pyleecan.Classes.MeshMat import MeshMat
from pyleecan.Classes.NodeMat import NodeMat
from pyleecan.Classes.CellMat import CellMat
from pyleecan.Classes.MeshSolution import MeshSolution
from pyleecan.Classes.SolutionMat import SolutionMat


def run_radial(
    self,
    axes_dict,
    Is_val=None,
    type_coord_sys=2,
    N_point_r=37,
    Kmesh_fineness=2,
    rotor_shift=8,
):
    """
    Solve the MagNetwork in the case of radial coordinate system type

    Parameters
    ----------
    self : MagNetwork
        A MagNetwork object
    axes_dict : dict
        Dict containing axes for MagNetwork calculation
    Is_val: ndarray
        Stator current matrix accounting for magnetic periodicities
    type_coord_sys : integer (Default = 2)
        Type of the coordinate system : 1 for Cartesian and 2 for Radial
    Kmesh_fineness : integer (Default = 2)
        Density of mesh cells inside the geometry
    rotor_shift : integer (Default = 8)
        Number of rotor mesh cells to be shifted with respect to the stator

    Returns
    -------
    Bx : nd-array
        Flux density in the radial direction
    By: nd-array
        Flux density in the theta direction
    Bx_airgap :
        Flux density of the airgap in the radial direction
    By_airgap :
        Flux density of the airgap in the theta direction
    """
    Machine = self.parent.machine
    # Geometry inputs
    la = Machine.rotor.L1  # Active length (m)
    Br = Machine.rotor.magnet.mat_type.mag.Brm20
    mu0 = np.pi * 4e-7  # Permeability of vacuum (H/m)

    # Definition of N_point_theta initial value
    if Machine.comp_periodicity_spatial()[1] == True:
        periodicity = Machine.comp_periodicity_spatial()[0]
    else:
        periodicity = Machine.comp_periodicity_spatial()[0] / 2

    angle_tp = (np.pi / periodicity) * (180 / np.pi)

    N_point_theta = Kmesh_fineness * round(0.5 * angle_tp) + 1

    # Definition of N_point_r
    # N_point_r = 1 + Kmesh_fineness * round(
    #     (Machine.stator.Rext - Machine.rotor.Rint) * 1000
    # )

    # Update of N_point_theta verifying condition 1
    N_point_theta = self.geometry_motor(N_point_theta, N_point_r, rotor_shift)[2]

    # Definition of the r-axis
    r = np.linspace(Machine.rotor.Rint, Machine.stator.Rext, N_point_r)

    # Definition of the theta axis
    theta = axes_dict["theta_primal"].get_values(is_smallestperiod=True)

    # Definition of the theta_dual and r_dual axes
    theta_dual = axes_dict["angle"].get_values(is_smallestperiod=True)
    r_dual = (r[1:] + r[:-1]) / 2
    # theta_dual = (theta[1:] + theta[:-1]) / 2

    # Definition of the boundary conditions
    BC = [
        "anti_periodic_condition",
        "homogeneous_Dirichlet_condition",
        "anti_periodic_condition",
        "homogeneous_Dirichlet_condition",
    ]

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

    # solving the model using the solver_linear_model
    (
        Phi,
        list_geometry,
        Num_unknowns,
        list_elem,
        permeability_cell,
        list_coord,
    ) = self.solver_linear_model(
        N_point_theta,
        N_point_r,
        theta,
        r,
        theta_dual,
        r_dual,
        rotor_shift,
        BC,
        self.geometry_motor,
        mu0,
        la,
        Br,
        type_coord_sys,
        JA=JA,
        JB=JB,
        JC=JC,
    )

    # Transfomration of radial coordinates to cartesian to plot the flux density contour
    x = (list_coord[:, 1] * np.cos(list_coord[:, 0])).reshape(N_point_r, N_point_theta)
    y = (list_coord[:, 1] * np.sin(list_coord[:, 0])).reshape(N_point_r, N_point_theta)

    list_cartesian_coord = np.zeros(list_coord.shape)
    list_cartesian_coord[:, 0] = x.flatten()
    list_cartesian_coord[:, 1] = y.flatten()

    # Plotting the flux density contour
    self.view_contour_flux(Phi, x, y, N_point_theta, N_point_r, list_geometry)

    # computing B radial
    Bx, By = self.compute_B(Phi, list_elem, list_coord, la, type_coord_sys)

    B = np.stack((Bx, By), axis=-1)

    temp = np.zeros((list_elem.shape[0], 3))
    temp[:, 0] = Bx.flatten()
    temp[:, 1] = By.flatten()
    B = temp

    temp = np.zeros((list_coord.shape[0], 3))
    temp[:, 0:2] = list_coord

    print("mesh saved", list_coord.shape, list_elem.shape)

    # Compute 2D curve of the airgap flux density
    Bx_airgap, By_airgap = self.comp_flux_airgap_local(
        # r,
        # theta,
        r_dual,
        theta_dual,
        Phi,
        list_elem,
        list_coord,
        la,
        Machine.comp_Rgap_mec(),
        type_coord_sys,
    )

    return Bx, By, Bx_airgap, By_airgap
