# -*- coding: utf-8 -*-


def init_mesh(
    self,
    N_point_theta,
    N_point_r,
    theta,
    r,
    boundary_condition,
):
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
