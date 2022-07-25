# -*- coding: utf-8 -*-


def init_mesh(
    self,
    N_point_theta,
    N_point_r,
    theta,
    r,
):
    list_coord = self.init_point(N_point_theta, N_point_r, theta, r)

    list_elem = self.init_cell(N_point_theta, N_point_r)

    return (
        list_coord,
        list_elem,
    )
