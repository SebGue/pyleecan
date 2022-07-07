# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:37:51 2022

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
from post_processing import add_BC_to_F, compute_B_square
from assembler import assembly, right_member_assembly, assembly_one_area
"""
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
    self,
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
    """
    Update permeability for non-linear model (saturated).

    Parameters
    ----------
    Phi : nd-array, size: k (floats)
        Flux with no boundary conditions.
    permeability_cell : nd-array, size: n (int)
        Permeability in each cell.
    list_geometry : nd-array, size: (foat)
        Materials in the cells.
    Num_Unknowns : nd-array, size: n (int)
        Numerotation in the linear system.
    list_elem : nd-array, size: m (integers)
        tab of elements
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.
    list_coord : nd-array of coordinate ( n x 2 )
        coordinate of each vertices
    la : float
        depth of the model.
    mu0 : float
        permeability of the vaccum.

    Returns
    -------
    permeability_cell : nd-array, size: n (int)
        Permeability in each cell.

    """
    # Compute staured mu
    B_sat = 1.99
    materials = np.array([5, 7])
    mu_r0 = 7500

    self.Phi = self.add_BC_to_F(
        self.Phi, self.Num_Unknowns, self.list_elem, self.BC_list
    )
    Bx, By = self.compute_B_square(self.Phi, self.list_elem, self.list_coord, self.la)
    B_norm = np.sqrt(Bx * Bx + By * By)

    H = B_norm / (self.mu0 * self.permeability_cell)
    mask = (self.list_geometry == materials[0]) + (self.list_geometry == materials[1])
    print(
        "H=",
        np.round(H[mask].max(), 5),
        "A/m, mu=",
        np.round(self.permeability_cell[mask].min(), 5),
        "H/m",
    )

    self.permeability_cell[mask] = 1 + (
        2 * B_sat / (np.pi * self.mu0 * H[mask])
    ) * np.arctan(np.pi * (mu_r0 - 1) * self.mu0 * H[mask] / (2 * B_sat))
    # print(permeability_cell[mask])
    return permeability_cell


def non_linear_model(self, size_x, size_y, x, y, pos, geometry, sol0, BC, la, mu0, Br):
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
    sol0 : nd-array
        first step.
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
    # Initialyze nature of elements
    h_x = (self.x.max() - self.x.min()) / (self.size_x - 1)
    h_y = (self.y.max() - self.y.min()) / (self.size_y - 1)
    list_geometry, permeability_materials = self.geometry(
        self.size_x, self.size_y, h_x, h_y, self.pos
    )

    non_linear_area = np.array([5, 7], dtype=np.int32)

    # Initialyze mesh
    list_coord = self.init_point(self.size_x, self.size_y, self.x, self.y)

    permeability_cell = self.init_permeabilty_cell(
        self.size_x, self.size_y, permeability_materials, mu0, list_geometry
    )

    list_elem = self.init_cell(self.size_x, self.size_y)

    BC_list, Periodic_point = self.init_mesh_BC(self.size_x, self.size_y, self.BC)

    Num_Unknowns = self.numeroting_unknows(list_elem, BC_list, Periodic_point)

    nn = Num_Unknowns.max() + 1

    reluc_list = self.init_reluc(list_elem, list_coord, self.mu0, self.la)

    t1 = time.perf_counter()
    print("Assembly geometry:", np.round(t1 - t0, 5), "secondes")
    # Assembly each area of the matrix
    t2 = time.perf_counter()

    t3 = time.perf_counter()

    # Assembly RHS
    E = self.right_member_assembly(
        list_geometry, Num_Unknowns, list_elem, list_coord, self.Br, self.mu0
    )

    t4 = time.perf_counter()
    print("Assembly vector:", np.round(t4 - t3, 5), "secondes")
    print("Total :", np.round(t4 - t2, 5))
    t12 = time.perf_counter()

    M_lin = csr_matrix((nn, nn), dtype=np.float64)
    M_non_lin = csr_matrix((nn, nn), dtype=np.float64)

    # Togather area of the matrix
    for i in range(0, permeability_materials.size):
        Matrix = self.assembly_one_area(
            i + 1,
            reluc_list,
            list_geometry,
            Num_Unknowns,
            list_elem,
            permeability_cell,
            BC_list,
        )
        print(nn, Matrix.shape)
        if (i + 1) in non_linear_area:
            M_non_lin += Matrix
        else:
            M_lin += Matrix

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
    criteria2 = eps + 1
    # Non linear loop

    while criteria2 > eps and i < 500:
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
            self.la,
            self.mu0,
        )
        # Update matrix
        M_non_lin = csr_matrix((nn, nn), dtype=np.float64)

        for j in non_linear_area:
            Matrix = self.assembly_one_area(
                j,
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

    F = self.add_BC_to_F(F, Num_Unknowns, list_elem, BC_list)

    return F, list_geometry, Num_Unknowns, list_elem, permeability_cell, list_coord
