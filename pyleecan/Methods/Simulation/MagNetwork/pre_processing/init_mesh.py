# -*- coding: utf-8 -*-


def init_mesh(
    self,
    N_point_theta,
    N_point_r,
    theta,
    r,
    boundary_condition,
):
    """Initialize the mesh of MagNetwork

    Parameters
    ----------
    self : MagNetwork
        the MagNetwork object to update
    N_point_theta : integer
        Number of discretization points in the theta-direction
    N_point_r : integer
        Number of discretization points in the r-direction
    theta: array
        Values of theta coordinates
    r: array
        Values of r coordinates
    boundary_condition : list of strings
        Boundary conditions of the simulation

    Returns
    -------
    list_coord : nd-array (size: N_point_r * N_point_theta * 2 of loats )
        list of point coordinates
    list_elem : nd-array, size: m (integers)
        tab of rectangular elements
    list_elem_permeability : nd-array (size: n of integers)
        Permeability in each cell inside the motor geometry
    boundary_condition_list : nd-array, size: m (int)
        Data strucuture
    Num_Unknowns :
        Numerotation in the linear system
    """
    list_coord = self.init_point(N_point_theta, N_point_r, theta, r)

    list_elem = self.init_cell(N_point_theta, N_point_r)

    list_elem_permability = self.geometry_motor(
        N_point_theta, N_point_r, self.rotor_shift
    )[0]

    list_boundary_condition, Periodic_point = self.init_mesh_BC(
        N_point_theta, N_point_r, boundary_condition
    )

    Num_Unknowns = self.numeroting_unknows(
        list_elem, list_boundary_condition, Periodic_point
    )

    return (
        list_coord,
        list_elem,
        list_elem_permability,
        list_boundary_condition,
        Num_Unknowns,
    )
