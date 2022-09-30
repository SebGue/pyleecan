# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 14:31:15 2022

@author: LAP02
"""


from pyleecan.Classes.MeshMat import MeshMat
from pyleecan.Classes.NodeMat import NodeMat
from pyleecan.Classes.CellMat import CellMat
from pyleecan.Classes.MeshSolution import MeshSolution
from pyleecan.Classes.SolutionMat import SolutionMat
import numpy as np


def cartesianmeshclass_pyleecan(self):

    Machine = self.parent.machine
    # Machine's characteristics
    la = Machine.rotor.L1  # Active length (m)

    # Br = 1.2  # PM remanent induction (residual induction) (T)
    Br = Machine.rotor.magnet.mat_type.mag.Brm20
    mu0 = np.pi * 4e-7  # Permeability of vacuum (H/m)
    mur1 = 1  # Relative permeability of air

    # mur2 = 7500  # Relative permeability of statot iron
    mur2 = Machine.stator.mat_type.mag.mur_lin

    # mur3 = 7500
    mur3 = Machine.rotor.mat_type.mag.mur_lin

    # Relative permeability of the winding
    if Machine.stator.winding.conductor.cond_mat.mag != None:
        mur_bob = Machine.stator.winding.conductor.cond_mat.mag.mur_lin
    else:
        mur_bob = 1

    # Relative permeabiltity of the PM
    mur_PM = Machine.rotor.magnet.mat_type.mag.mur_lin

    list_materials = ["bob1", "bob2", "bob3", "air", "iron1", "PM", "iron3"]

    # permeabilty_materials = np.array([1, 1, 1, 1, 7500, 1, 7500])
    permeabilty_materials = np.array(
        [mur_bob, mur_bob, mur_bob, mur1, mur2, mur_PM, mur3]
    )

    x_min = 0
    x_max = (
        np.pi * Machine.stator.Rint / Machine.rotor.get_pole_pair_number()
    )  # equal to tp, in meters

    y_min = 0
    y_max = Machine.stator.Rext  # in meters

    mul = 2
    size_x = int(x_max * 1000 * mul) + 1
    size_y = int(y_max * 1000 * mul) + 1

    x = np.linspace(x_min, x_max, size_x)
    y = np.linspace(y_min, y_max, size_y)

    x_dual = (x[1:] + x[: size_x - 1]) / 2
    y_dual = (y[1:] + y[: size_y - 1]) / 2

    BC = ["AP", "HD", "AP", "HD"]
    pos = 0
    mode = "cartesian"

    (
        F,
        list_geometry,
        Num_unknowns,
        list_elem,
        permeability_cell,
        list_coord,
    ) = self.solver_linear_model(
        size_x,
        size_y,
        x,
        y,
        x_dual,
        y_dual,
        pos,
        BC,
        self.geometry_motor,
        mu0,
        la,
        Br,
        mode,
        JA=None,
        JB=None,
        JC=None,
    )

    Bx, By = self.compute_B_square(F, list_elem, list_coord, la)
    B = np.zeros((Bx.size, 2))
    B[:, 0] = Bx
    B[:, 1] = By
    Norm_B = np.sqrt(Bx ** 2 + By ** 2)

    print("Solve RN done.")
    mesh = MeshMat(dimension=3)
    mesh.node = NodeMat()
    print("Add points in mesh")
    for i in range(list_coord.shape[0]):
        mesh.node.add_node([list_coord[i, 0], list_coord[i, 1], 0])
    print("Done \nAdd elements in mesh")

    mesh.cell["quad"] = CellMat(nb_node_per_cell=4)
    for i in range(list_elem.shape[0]):
        mesh.add_cell(list_elem[i, :], "quad")

    print("Done \nAdd material for ech elementss")

    MSol = MeshSolution(mesh=[mesh])

    for i in range(list_elem.shape[0]):
        MSol.group = {list_materials[list_geometry[i] - 1]: list_elem[i, :]}

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

    # plot B.

    # field_vec=B[np.newaxis]
    # print(field_vec.shape)# (1, 12240, 2)
    # my_vec_solution = SolutionMat(
    #     label="my_vector",
    #     type_cell="cell",
    #     field=field_vec,
    #     indice=np.arange(list_elem.shape[0]),
    #     axis_name=["time", "indice", "component"],
    #     axis_size = [1,B.shape[0],B.shape[1]]
    # )
    # MSol.solution.append(my_vec_solution)
    # MSol.plot_glyph(label="my_vector", is_point_arrow=True, factor=1/10)
