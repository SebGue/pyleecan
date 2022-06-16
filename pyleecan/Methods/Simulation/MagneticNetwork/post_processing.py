# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:11:39 2022

@author: LAP02
"""

import numpy as np


def compute_B_square(F, list_elem,list_coord,la):
    """
    Compute the magnetic field from the flux F and the mesh.

    Parameters
    ----------
    F : nd-array, size: n (float)
        Flux
    list_elem : nd-array, size: m (integers)
        tab of elements
    list_coord : nd-array, size: nx2 (flat)
        coordinate of each point.
    la : float
        depth of the model.

    Returns
    -------
    Bx : nd-array, size: m (flat)
        x component of B
    By : nd-array, size: m (float)
        y component of B

    """
    
    h_x=np.linalg.norm(list_coord[list_elem[:, 0]]-list_coord[list_elem[:, 1]],axis=1,ord=2)
    h_y=np.linalg.norm(list_coord[list_elem[:, 1]]-list_coord[list_elem[:, 2]],axis=1,ord=2)
    
    
    # Compute B
    Bx1 = (F[list_elem[:, 2]]-F[list_elem[:, 1]])/(h_y*la)
    Bx2 = (F[list_elem[:, 3]]-F[list_elem[:, 0]])/(h_y*la)
    
    By1 = (F[list_elem[:, 0]]-F[list_elem[:, 1]])/(h_x*la)
    By2 = (F[list_elem[:, 3]]-F[list_elem[:, 2]])/(h_x*la)


    Bx = ((Bx1+Bx2)/2)
    By = ((By1+By2)/2)

    return Bx, By


def add_BC_to_F(F,Num_Unknowns, list_elem,BC_list):
    """
    Add boundary condition to the the flux.

    Parameters
    ----------
    F : nd-array, size: Num_Unknows.max()+1 (float)
        Flux
    Num_Unknowns : nd-array, size: n (integers)
        list of unknowns in the linear system .
    list_elem : nd-array, size: m (integers)
        tab of elements
    BC_list : nd-array, size: n (integers) 
        data-strurture for BC evaluation.

    Returns
    -------
    F_full : size: n (float)
        The flux on each point of the mesh.

    """
    #post processing procedure to add boundary condition 
    
    F_full= np.zeros(Num_Unknowns.size)
    mask=Num_Unknowns!=-1
    
    F_full[mask]=F[Num_Unknowns[mask]]

    if 3 in BC_list:
        mask_AP=BC_list == 3
        F_full[mask_AP] = -F_full[mask_AP]

    return F_full