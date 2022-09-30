# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:37:51 2022

@author: LAP02
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
#Choose beetween chol-mod and scipy sparse linear solver
Have_cholmod=False
try:
    from sksparse.cholmod import cholesky,analyze
    Have_cholmod=True
except:
    from scipy.sparse.linalg import spsolve



def mu_non_linear(self,Phi, permeability_cell, list_geometry,Num_Unknowns,list_elem,BC_list,list_coord,la,mu0):
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
    materials = np.array([5,7])
    mu_r0 = 7500
    
    
    Phi=self.add_BC_to_Phi(Phi, Num_Unknowns,list_elem,BC_list)
    Bx, By = self.compute_B_square(Phi, list_elem,list_coord,la)
    B_norm = np.sqrt(Bx*Bx+By*By)

    H = B_norm/(mu0*permeability_cell)
    mask = (list_geometry == materials[0])+(list_geometry == materials[1])
    print("H=", np.round(H[mask].max(),5), "A/m, mu=", np.round(permeability_cell[mask].min(),5),"H/m")

    permeability_cell[mask] = 1+(2*B_sat/(np.pi*mu0*H[mask])) * \
        np.arctan(np.pi*(mu_r0-1)*mu0*H[mask]/(2*B_sat))
    # print(permeability_cell[mask])
    return permeability_cell


def non_linear_model(self,
    N_point_theta,
    N_point_r,
    theta,
    r,
    x_dual,
    y_dual,
    boundary_condition,
    geometry,
    material_dict,
    la,
    Br,
    type_coord_sys,
    JA=None,
    JB=None,
    JC=None,
    sol0=None):
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
    non_linear_area=np.array([5,7],dtype=np.int32)
    
    # Initialyze mesh
    (
        list_coord,
        list_elem,
        list_elem_permability,
        list_boundary_condition,
        Num_Unknowns,
    ) = self.init_mesh(
        N_point_theta,
        N_point_r,
        theta,
        r,
        boundary_condition,
    )
    
    nn=Num_Unknowns.max()+1

    mu0 = material_dict["vacuum"]
    reluc_list = self.init_reluc(list_elem, list_coord, mu0, la, type_coord_sys)
    
    t1 = time.perf_counter()
    print("Assembly geometry+mesh:", np.round(t1-t0, 5), "secondes")
    # Assembly each area of the matrix
    t2 = time.perf_counter()


    t3 = time.perf_counter()


    # Assembly RHS
    RHS = self.right_member_assembly(
        list_elem_permability,
        Num_Unknowns,
        list_elem,
        list_coord,
        reluc_list,
        Br,
        material_dict,
        None,
        la,
        type_coord_sys,
        N_point_theta,
        JA=JA,
        JB=JB,
        JC=JC,
    )

    
    
    t4 = time.perf_counter()
    print("Assembly vector:", np.round(t4-t3, 5), "secondes")
    print("Total :", np.round(t4-t2, 5))
    t12 = time.perf_counter()
    
    M_lin= csr_matrix((nn,nn), dtype=np.float64)
    M_non_lin= csr_matrix((nn,nn), dtype=np.float64)
    
    # Togather area of model with tghe matrix
    for i in range(0,list_elem_permability.size):
        Matrix=self.assembly_one_area(i+1, reluc_list,list_boundary_condition, Num_Unknowns, list_elem, permeability_cell,list_boundary_condition)
        print(nn,Matrix.shape)
        if (i+1) in non_linear_area:
            M_non_lin += Matrix
        else:
            M_lin += Matrix
        

    M_csr=M_non_lin+M_lin

    # Initialize non linear loop
    if sol0==None:
        if Have_cholmod:
            M_csc=M_csr.tocsc()
            #Symbolic Cholesky factorisation
            factor=analyze(M_csc,mode="simplicial",ordering_method="nesdis")
            #Cholesky factorisation
            factor.cholesky_inplace(M_csc)
            Phi = factor(RHS)
        else:
            Phi = spsolve(M_csr, RHS, permc_spec="MMD_AT_PLUS_A", use_umfpack=True)
    else:
        Phi=sol0
    
    t1 = time.perf_counter()
    old_permaebility = permeability_cell.copy()
    phi_old = Phi.copy()
    i = 0
    eps = 10**-4
    criteria = eps+1
    criteria2 = eps+1
    # Non linear loop
    
    while criteria2 > eps and i < 500:
        print("loop n°", i+1)
        # Update permeability
        permeability_cell = mu_non_linear(Phi, permeability_cell, material_dict,Num_Unknowns,list_elem,list_boundary_condition,list_coord,la,mu0)
        # Update matrix
        M_non_lin= csr_matrix((nn,nn), dtype=np.float64)
        
        for j in material_dict:
            Matrix=self.assembly_one_area(j, reluc_list,material_dict, Num_Unknowns, list_elem, permeability_cell,list_boundary_condition)
            M_non_lin += Matrix
            
            
        M_csr=M_non_lin+M_lin
        
        if Have_cholmod:
           # Compute cholesky decomposition, solve the system, compute solution
           factor.cholesky_inplace( M_csr.tocsc())
           Phi = factor(RHS)
        else:
           Phi = spsolve(M_csr, RHS, permc_spec="MMD_AT_PLUS_A", use_umfpack=True) 
        
       
        # Stoping criterion
        criteria = np.linalg.norm(Phi-phi_old, ord=2)/np.linalg.norm(Phi, ord=2)
        criteria2= np.linalg.norm(old_permaebility-permeability_cell, ord=2)/np.linalg.norm(permeability_cell, ord=2)
        # Updtae old solution
        phi_old = Phi.copy()
        old_permaebility = permeability_cell.copy()
        i += 1
        print("eps flux:", np.round(criteria,7),"eps permeability", np.round(criteria2,7) ,"\n")
    
    t2 = time.perf_counter()
    print("Time to solve non-linear:", np.round(t2-t1, 5),
          "secondes, res:", np.linalg.norm( M_csr@Phi-RHS, ord=2))
    
    Phi=self.add_BC_to_Phi(Phi, Num_Unknowns,list_elem,list_boundary_condition)
    
    
    return (
        Phi,
        Num_Unknowns,
        list_elem,
        list_elem_permability,
        list_coord,
    )
