# -*- coding: utf-8 -*-


def init_mesh(
    self,
    N_point_theta,
    N_point_r,
    theta,
    r,
    permeability_materials,
    mu0,
    list_elem_materials,
    boundary_condition,
):
    list_coord = self.init_point(N_point_theta, N_point_r, theta, r)

    list_elem_permability = self.init_permeabilty_cell(
        N_point_theta, N_point_r, permeability_materials, list_elem_materials
    )

    list_elem = self.init_cell(N_point_theta, N_point_r)

    list_boundary_condition, Periodic_point = self.init_mesh_BC(
        N_point_theta, N_point_r, boundary_condition
    )

    Num_Unknowns = self.numeroting_unknows(
        list_elem, list_boundary_condition, Periodic_point
    )

    self.save_mesh(
        list_elem_materials, Num_Unknowns, list_elem, theta, r, list_boundary_condition
    )

    return (
        list_coord,
        list_elem_permability,
        list_elem,
        list_boundary_condition,
        Periodic_point,
        Num_Unknowns,
    )
