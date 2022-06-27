# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:10:31 2022

@author: LAP02
"""

"""
from pre_processing import (
    save_mesh,
    init_reluc,
    init_point,
    init_cell,
    init_permeabilty_cell,
    init_permeabilty_cell,
    init_mesh_BC,
    numeroting_unknows,
)
from post_processing import add_BC_to_F
from assembler import assembly, right_member_assembly
"""
import time

# import geometry_linear_motor

from threadpoolctl import threadpool_info
import os

# from pyleecan.Classes.MagneticNetwork import MagneticNetwork

threads = "1"
os.environ["OMP_NUM_THREADS"] = threads  # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = threads  # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = threads  # export MKL_NUM_THREADS=6
# export VECLIB_MAXIMUM_THREADS=4
os.environ["VECLIB_MAXIMUM_THREADS"] = threads
os.environ["NUMEXPR_NUM_THREADS"] = threads  # export NUMEXPR_NUM_THREADS=6
# export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=2"
os.environ["MKL_DOMAIN_NUM_THREADS"] = "MKL_DOMAIN_ALL=1"


os.environ["MKL_DYNAMIC"] = "FALSE"
os.environ["OMP_DYNAMIC"] = "FALSE"

import numpy as np

# Choose beetween chol-mod and scipy sparse linear solver
Have_cholmod = False
try:
    from sksparse.cholmod import cholesky

    Have_cholmod = True
except:
    from scipy.sparse.linalg import spsolve


def solver_linear_model(
    self, size_x, size_y, x, y, x_dual, y_dual, pos, BC, geometry, mu0, la, Br, mode
):
    """

    Parameters
    ----------
    size_x : integer
        Size of x.
    size_y : integer
        Size of y.
    x : nd-array, size: size_x (float)
        x coordinate.
    y : nd-array, size: size_x (float)
        y coordinate.
    pos : integers
        position of the rotor.
    geometry : object function
        function to initialyze model.
    BC : list of string (condition)
        List of boundary conditions.
    la : float
        depth of the model.
    mu0 : float
        permeability of the vaccum.
    Br : float
        Caracterization of the permanent magnet.

    Returns
    -------
    F : size: n (float)
        The flux on each point of the mesh.
    list_geomerty: nd-array, size: m (integers)
        contains the material of each cells (elements)
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    permeability_cell : nd-array, size: m (float)
        the permeability values of each cell
    list_coord : nd-array, size: size_x*size_y x 2 (float)
        list of coordinate.

    """

    t0 = time.perf_counter()
    import matplotlib.pyplot as plt

    # Initialyze nature of elements
    list_geometry, permeability_materials = self.geometry_linear_motor(
        size_x, size_y, pos
    )

    list_coord = self.init_point(size_x, size_y, x, y)
    print("Permeability", permeability_materials)
    permeability_cell = self.init_permeabilty_cell(
        size_x, size_y, permeability_materials, list_geometry
    )

    list_elem = self.init_cell(size_x, size_y)

    BC_list, Periodic_point = self.init_mesh_BC(size_x, size_y, BC)

    Num_Unknowns = self.numeroting_unknows(list_elem, BC_list, Periodic_point)

    t1 = time.perf_counter()
    print("Assembly geometry:", np.round(t1 - t0, 5), "seconds")
    self.save_mesh(list_geometry, Num_Unknowns, list_elem, x, y, BC_list)
    t2 = time.perf_counter()
    print("Save mesh:", np.round(t2 - t1, 5), "secondes")

    # Assembly all matrice
    reluc_list = self.init_reluc(list_elem, list_coord, mu0, la, mode)
    print(reluc_list)
    M_csr = self.assembly(
        reluc_list, Num_Unknowns, list_elem, permeability_cell, BC_list
    )

    t3 = time.perf_counter()
    print("Assembly matrix", np.round(t3 - t2, 5), "secondes")

    # Assembly RHS
    E = self.right_member_assembly(
        list_geometry, Num_Unknowns, list_elem, list_coord, Br, mu0, mode
    )

    t4 = time.perf_counter()

    print("Assembly vector:", np.round(t4 - t3, 5), "secondes")
    print("Total :", np.round(t4 - t2, 5))

    # Compute Solution
    t3 = time.perf_counter()
    if Have_cholmod:

        # Compute the solution
        factor = cholesky(M_csr.tocsc())
        F = factor(E)
        t4 = time.perf_counter()
        print(
            "Time to solve (CholMod):",
            np.round(t4 - t3, 5),
            "secondes, res:",
            np.linalg.norm(M_csr @ F - E, ord=2),
        )
    else:
        F = spsolve(M_csr, E, permc_spec="MMD_AT_PLUS_A", use_umfpack=True)
        t4 = time.perf_counter()
        print(
            "Time to solve (direct,UMF):",
            np.round(t4 - t3, 5),
            "secondes, res:",
            np.linalg.norm(M_csr @ F - E, ord=2),
        )

    F = self.add_BC_to_F(F, Num_Unknowns, list_elem, BC_list)

    return F, list_geometry, Num_Unknowns, list_elem, permeability_cell, list_coord
