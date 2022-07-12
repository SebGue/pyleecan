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

from pyleecan.Classes.MeshMat import MeshMat
from pyleecan.Classes.NodeMat import NodeMat
from pyleecan.Classes.CellMat import CellMat
from pyleecan.Classes.MeshSolution import MeshSolution
from pyleecan.Classes.SolutionMat import SolutionMat


def run_radial(self, axes_dict, Is_val=None):

    Machine = self.parent.machine
    la = Machine.rotor.L1  # Active length (m)

    # Br = 1.2  # PM remanent induction (residual induction) (T)
    Br = Machine.rotor.magnet.mat_type.mag.Brm20
    mu0 = np.pi * 4e-7  # Permeability of vacuum (H/m)

    pos = 2

    x_min = 0

    # Theta angle represented by x_max
    x_max = np.pi / Machine.comp_periodicity_spatial()[0]

    # Definition of the y-axis
    y_min = 0
    y_max = Machine.stator.Rext - Machine.rotor.Rint  # in meters

    # Density of the mesh
    density = 1

    # Size of the mesh according to r
    size_r = (int)(density * y_max * 1000 + 1)  # step of discretization of r

    # Size of mesh according to theta
    size_theta = (int)(density * x_max * 1000 + 1)  # step of discretization of theta

    # Definition of the r-axis
    r = np.linspace(Machine.rotor.Rint, Machine.stator.Rext, size_r)

    # Definition of the theta-axis
    # Add one extra point so that dual mesh has the correct dimension
    theta = np.linspace(
        axes_dict["angle"].initial,
        axes_dict["angle"].final,
        axes_dict["angle"].number + 1,
        endpoint=False,
    )

    # angle = np.pi / Machine.comp_periodicity_spatial()[0]
    # theta = np.linspace(0, angle, size_theta)

    # Definition of the r_dual of the dual mesh
    r_dual = (r[1:] + r[:-1]) / 2

    # Definition of the theta_dual of the dual mesh
    theta_dual = (theta[1:] + theta[-1]) / 2

    # Definition of the boundary conditions list
    BC = ["AP", "HD", "AP", "HD"]

    # Definition of the mode, polar in the case of radial geometry
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

    list_cart = np.zeros((list_coord.shape[0], 2))

    list_cart[:, 0] = x.flatten()
    list_cart[:, 1] = y.flatten()

    self.view_contour_flux(F, x, y, x.shape[1], x.shape[0], list_geometry)

    # print(self.view_contour_flux(F, x, y, x.shape[1], x.shape[0], list_geometry))

    Bx, By = self.compute_B_square(F, list_elem, list_coord, la)

    B = np.stack((Bx, By), axis=-1)

    temp = np.zeros((list_elem.shape[0], 3))
    temp[:, 0] = Bx.flatten()
    temp[:, 1] = By.flatten()
    B = temp

    temp = np.zeros((list_coord.shape[0], 3))
    temp[:, 0:2] = list_coord

    # Launch mesh.io
    # list_coord = temp
    # points = list_coord
    # cells = [
    #     ("quad", list_elem),
    # ]
    # mesh = meshio.Mesh(
    #     points,
    #     cells,
    #     # Optionally provide extra data on points, cells, etc.
    #     point_data={"Flux": F},
    #     # Each item in cell data must match the cells array
    #     cell_data={"B": [B], "Materials": [list_geometry]},
    # )
    # mesh.write(
    #     "mymesh.xdmf",  # str, os.PathLike, or buffer/open file
    #     # file_format="vtk",  # optional if first argument is a path; inferred from extension
    # )

    # Alternative with the same options
    # meshio.write_points_cells("mymesh.vtu", points, cells)

    print("mesh saved", list_coord.shape, list_elem.shape)

    # Compute 2D curve of the airgap flux density
    Bx_airgap, By_airgap = self.comp_flux_airgap_local(
        r, theta, F, list_elem, list_coord, la, Machine.comp_Rgap_mec()
    )
    # Add my mesh to pyleecan
    print("Solve RN done.")
    mesh = MeshMat(dimension=3)
    mesh.node = NodeMat()
    print("Add points in mesh")
    for i in range(list_cart.shape[0]):
        mesh.node.add_node([list_cart[i, 0], list_cart[i, 1], 0])
    print("Done \nAdd elements in mesh")

    mesh.cell["quad"] = CellMat(nb_node_per_cell=4)
    for i in range(list_elem.shape[0]):
        mesh.add_cell(list_elem[i, :], "quad")

    MSol = MeshSolution(mesh=[mesh])

    # print("Done \nAdd material for ech elementss")
    # for i in range(list_elem.shape[0]):
    #     MSol.group = {list_materials[list_geometry[i] - 1]: list_elem[i, :]}

    print("Done")
    # plot the mesh
    MSol.plot_mesh()

    # plot the flux
    field = F[np.newaxis]
    print(field.shape)
    my_solution = SolutionMat(
        label="Flux (Weber)",
        type_cell="point",
        field=field,
        indice=np.arange(list_coord.shape[0]),
        axis_name=["time", "indice"],
        axis_size=[1, list_coord.shape[0]],
    )
    MSol.solution.append(my_solution)
    MSol.plot_contour()

    return Bx, By, Bx_airgap, By_airgap
