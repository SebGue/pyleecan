# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:10:31 2022
@author: LAP02
"""
import time

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
import matplotlib.pyplot as plt

# Choose beetween chol-mod and scipy sparse linear solver
Have_cholmod = False
try:
    from sksparse.cholmod import cholesky

    Have_cholmod = True
except:
    from scipy.sparse.linalg import spsolve


def solver_linear_model(
    self,
    N_point_theta,
    N_point_r,
    theta,
    r,
    x_dual,
    y_dual,
    rotor_shift,
    BC,
    geometry,
    mu0,
    la,
    Br,
    mode,
    JA=None,
    JB=None,
    JC=None,
):
    """
    Parameters
    ----------
    N_point_theta : integer
        Size of theta.
    N_point_r : integer
        Size of r.
    theta : nd-array, size: N_point_theta (float)
        theta coordinate.
    r : nd-array, size: N_point_theta (float)
        r coordinate.
    rotor_shift : integers
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
    Phi : size: N_point_theta (float)
        The flux on each point of the mesh.
    list_geomerty: nd-array, size: m (integers)
        contains the material of each cells (elements)
    Num_Unknowns : nd-array, size: N_point_theta (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    list_elem_permability : nd-array, size: m (float)
        the permeability values of each cell
    list_coord : nd-array, size: N_point_theta*N_point_r theta 2 (float)
        list of coordinate.
    """

    t0 = time.perf_counter()

    # initialize the cells of materials and their permeability
    list_elem_materials, permeability_materials = self.geometry_motor(
        N_point_theta, N_point_r, rotor_shift
    )

    # initialize the list_coord which contains the grid of points
    list_coord = self.init_point(N_point_theta, N_point_r, theta, r)
    # print("Permeability", permeability_materials)

    list_elem_permability = self.init_permeabilty_cell(
        N_point_theta, N_point_r, permeability_materials, list_elem_materials
    )

    list_elem = self.init_cell(N_point_theta, N_point_r)

    list_boundary_condition, Periodic_point = self.init_mesh_BC(
        N_point_theta, N_point_r, BC
    )

    Num_Unknowns = self.numeroting_unknows(
        list_elem, list_boundary_condition, Periodic_point
    )

    # Mesuring the performance time
    t1 = time.perf_counter()

    # print("Assembly geometry:", np.round(t1 - t0, 5), "seconds")
    self.save_mesh(
        list_elem_materials, Num_Unknowns, list_elem, theta, r, list_boundary_condition
    )

    t2 = time.perf_counter()
    # Saving mesh time
    print("Save mesh:", np.round(t2 - t1, 5), "secondes")

    # Assembly all matrice
    reluc_list = self.init_reluc(list_elem, list_coord, mu0, la, mode)
    # print(reluc_list)

    M_csr = self.assembly(
        reluc_list,
        Num_Unknowns,
        list_elem,
        list_elem_permability,
        list_boundary_condition,
    )

    t3 = time.perf_counter()
    print("Assembly matrix", np.round(t3 - t2, 5), "secondes")

    # Assembly RHS containing the sources
    RHS = self.right_member_assembly(
        list_elem_materials,
        Num_Unknowns,
        list_elem,
        list_coord,
        reluc_list,
        permeability_materials,
        Br,
        mu0,
        la,
        mode,
        JA=JA,
        JB=JB,
        JC=JC,
    )

    t4 = time.perf_counter()

    print("Assembly vector:", np.round(t4 - t3, 5), "secondes")
    print("Total :", np.round(t4 - t2, 5))

    # Compute Solution
    t3 = time.perf_counter()
    if Have_cholmod:

        # Compute the solution
        factor = cholesky(M_csr.tocsc())
        Phi = factor(RHS)
        t4 = time.perf_counter()
        print(
            "Time to solve (CholMod):",
            np.round(t4 - t3, 5),
            "secondes, res:",
            np.linalg.norm(M_csr @ Phi - RHS, ord=2),
        )
    else:
        Phi = spsolve(M_csr, RHS, permc_spec="MMD_AT_PLUS_A", use_umfpack=True)
        t4 = time.perf_counter()
        print(
            "Time to solve (direct,UMF):",
            np.round(t4 - t3, 5),
            "secondes, res:",
            np.linalg.norm(M_csr @ Phi - RHS, ord=2),
        )

    Phi = self.add_BC_to_Phi(Phi, Num_Unknowns, list_elem, list_boundary_condition)

    return (
        Phi,
        list_elem_materials,
        Num_Unknowns,
        list_elem,
        list_elem_permability,
        list_coord,
    )
