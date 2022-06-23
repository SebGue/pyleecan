# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:37:51 2022

@author: LAP02
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
from post_processing import add_BC_to_F, compute_B_square
from assembler import assembly, right_member_assembly, assembly_one_area
import time

from threadpoolctl import threadpool_info
import os

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
from scipy.sparse import csr_matrix

# Choose beetween chol-mod and scipy sparse linear solver
Have_cholmod = False
try:
    from sksparse.cholmod import cholesky, analyze

    Have_cholmod = True
except:
    from scipy.sparse.linalg import spsolve


def mu_non_linear(
    Phi,
    permeability_cell,
    list_geometry,
    Num_Unknowns,
    list_elem,
    BC_list,
    list_coord,
    la,
    mu0,
):
    # Compute staured mu
    B_sat = 1.8
    materials = 2
    mu_r0 = 7500
    Phi = add_BC_to_F(Phi, Num_Unknowns, list_elem, BC_list)
    Bx, By = compute_B_square(Phi, list_elem, list_coord, la)
    B_norm = np.sqrt(Bx * Bx + By * By)

    H = B_norm / (mu0 * permeability_cell)
    mask = list_geometry == materials
    print(
        "H=",
        np.round(H[mask].max(), 5),
        "A/m, mu=",
        np.round(permeability_cell[mask].min(), 5),
        "H/m",
    )

    permeability_cell[mask] = 1 + (2 * B_sat / (np.pi * mu0 * H[mask])) * np.arctan(
        np.pi * (mu_r0 - 1) * mu0 * H[mask] / (2 * B_sat)
    )
    # print(permeability_cell[mask])
    return permeability_cell


def non_linear_model(size_x, size_y, x, y, pos, geometry, sol0, BC, la, mu0, Br):
    t0 = time.perf_counter()
    # Initialyze nature of elements
    h_x = (x.max() - x.min()) / (size_x - 1)
    h_y = (y.max() - y.min()) / (size_y - 1)
    list_geometry, permeability_materials = geometry(size_x, size_y, h_x, h_y, pos)
    non_linear_area = np.array([4, 6])

    # Initialyze mesh
    list_coord = init_point(size_x, size_y, x, y)

    permeability_cell = init_permeabilty_cell(
        size_x, size_y, permeability_materials, list_geometry
    )

    list_elem = init_cell(size_x, size_y, x, y)

    BC_list, Periodic_point = init_mesh_BC(size_x, size_y, x, y, BC)

    Num_Unknowns = numeroting_unknows(list_elem, BC_list, Periodic_point)
    nn = (Num_Unknowns >= 0).sum()

    reluc_list = init_reluc(list_elem, list_coord, mu0, la)

    t1 = time.perf_counter()
    print("Assembly geometry:", np.round(t1 - t0, 5), "secondes")
    # Assembly each area of the matrix
    t2 = time.perf_counter()
    Matrices = []

    for i in range(1, permeability_materials.size + 1):
        M_csr = assembly_one_area(
            i,
            reluc_list,
            list_geometry,
            Num_Unknowns,
            list_elem,
            permeability_cell,
            BC_list,
        )
        Matrices.append(M_csr)

    t3 = time.perf_counter()
    print("Assembly matrix", np.round(t3 - t2, 5), "secondes")

    # Assembly RHS
    E = right_member_assembly(
        list_geometry, Num_Unknowns, list_elem, list_coord, Br, mu0
    )

    t4 = time.perf_counter()
    print("Assembly vector:", np.round(t4 - t3, 5), "secondes")
    print("Total :", np.round(t4 - t2, 5))
    t12 = time.perf_counter()

    M_lin = csr_matrix((nn, nn), dtype=np.float64)
    M_non_lin = csr_matrix((nn, nn), dtype=np.float64)

    # Togather area of the matrix
    for i in range(0, permeability_materials.size):
        if i in non_linear_area:
            M_non_lin += Matrices[i]
        else:
            M_lin += Matrices[i]

    M_csr = M_non_lin + M_lin

    # Initialize non linear loop

    if Have_cholmod:
        M_csc = M_csr.tocsc()
        # Symbolic Cholesky factorisation
        factor = analyze(M_csc, mode="simplicial", ordering_method="nesdis")
        # Cholesky factorisation
        factor.cholesky_inplace(M_csc)
        F = factor(E)
    else:
        F = spsolve(M_csr, E, permc_spec="MMD_AT_PLUS_A", use_umfpack=True)

    t1 = time.perf_counter()
    old_permaebility = permeability_cell.copy()
    phi_old = F.copy()
    i = 0
    eps = 10 ** -4
    criteria = eps + 1
    # Non linear loop

    while criteria > eps and i < 500:
        print("loop n°", i + 1)
        # Update permeability
        permeability_cell = mu_non_linear(
            F,
            permeability_cell,
            list_geometry,
            Num_Unknowns,
            list_elem,
            BC_list,
            list_coord,
            la,
            mu0,
        )
        # Update matrix
        M_non_lin = csr_matrix((nn, nn), dtype=np.float64)
        for i in non_linear_area:
            Matrix = assembly_one_area(
                i,
                reluc_list,
                list_geometry,
                Num_Unknowns,
                list_elem,
                permeability_cell,
                BC_list,
            )
            M_non_lin += Matrix

        M_csr = M_non_lin + M_lin

        if Have_cholmod:
            # Compute cholesky decomposition, solve the system, compute solution
            factor.cholesky_inplace(M_csr.tocsc())
            F = factor(E)
        else:
            F = spsolve(M_csr, E, permc_spec="MMD_AT_PLUS_A", use_umfpack=True)

        # Stoping criterion
        criteria = np.linalg.norm(F - phi_old, ord=2) / np.linalg.norm(F, ord=2)
        criteria2 = np.linalg.norm(
            old_permaebility - permeability_cell, ord=2
        ) / np.linalg.norm(permeability_cell, ord=2)
        # Updtae old solution
        phi_old = F.copy()
        old_permaebility = permeability_cell.copy()
        i += 1
        print(
            "eps flux:",
            np.round(criteria, 7),
            "eps permeability",
            np.round(criteria2, 9),
            "\n",
        )

    t2 = time.perf_counter()
    print(
        "Time to solve non-linear:",
        np.round(t2 - t1, 5),
        "secondes, res:",
        np.linalg.norm(M_csr @ F - E, ord=2),
    )

    F = add_BC_to_F(F, size_x, size_y, Num_Unknowns, list_elem, BC_list)

    return F, list_geometry, Num_Unknowns, list_elem, permeability_cell, list_coord
