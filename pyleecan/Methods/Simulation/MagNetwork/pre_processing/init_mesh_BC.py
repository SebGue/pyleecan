# -*- coding: utf-8 -*-

import numpy as np


def init_mesh_BC(self, N_point_theta, N_point_r, boundary_condition):
    """
    Init data structure to compute boundary condition in the sytem.

    Parameters
    ----------
    N_point_theta : integer
        Number of horizontal points.
    N_point_r : integer
        Number of vertical point.
    boundary_condition : list of string (condition)
        List of boundary conditions.

    Returns
    -------
    boundary_condition_list : nd-array, size: m (int)
        Data strucuture.
    Periodic_point : nd-array, size: ? (int)
        Periodic point connection.

    """
    # boundary_condition: Boundary Condition, ["condi0","cond1","cond2","cond3"]
    # -> "homogeneous_Dirichlet_condition": homogeneous Dicrichlet conditions (=0)
    # -> "periodic_condition" : Periodic contitions
    # -> "anti_periodic_condition": Anti Periodic Condition
    # -> "non_attribuate_condition" :
    # boundary_condition_list -> Boudary condition

    N_total_point = N_point_theta * N_point_r
    ###########################################

    boundary_condition_list = np.zeros(N_total_point, dtype=np.uint8)
    Periodic_point = np.zeros(0, dtype=np.uint16)

    if np.all(
        boundary_condition
        == [
            "periodic_condition",
            "homogeneous_Dirichlet_condition",
            "periodic_condition",
            "homogeneous_Dirichlet_condition",
        ]
    ):

        # Initialyze the boundary_condition list

        # Vertical boundary_condition
        boundary_condition_list[::N_point_theta] = 2
        boundary_condition_list[N_point_theta - 1 :: N_point_theta] = 2

        # Horizontal boundary_condition
        boundary_condition_list[:N_point_theta] = 1
        boundary_condition_list[N_point_theta * (N_point_r - 1) : N_total_point] = 1

        # Connexion of periodic points
        Periodic_point = -np.ones((N_point_r - 2, 2), dtype=np.int32)
        Periodic_point[:, 0] = (
            np.arange(0, N_point_theta * (N_point_r - 2), N_point_theta) + N_point_theta
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + N_point_theta - 1

    elif np.all(
        boundary_condition
        == [
            "anti_periodic_condition",
            "homogeneous_Dirichlet_condition",
            "anti_periodic_condition",
            "homogeneous_Dirichlet_condition",
        ]
    ):
        # Initialyze the boundary_condition list
        # Vertical boundary_condition
        boundary_condition_list[::N_point_theta] = 2
        boundary_condition_list[N_point_theta - 1 :: N_point_theta] = 3
        # Horizontal boundary_condition
        boundary_condition_list[:N_point_theta] = 1
        boundary_condition_list[N_point_theta * (N_point_r - 1) : N_total_point] = 1

        Periodic_point = -np.ones((N_point_r - 2, 2), dtype=np.int32)
        Periodic_point[:, 0] = (
            np.arange(0, N_point_theta * (N_point_r - 2), N_point_theta) + N_point_theta
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + N_point_theta - 1

    elif np.all(
        boundary_condition
        == [
            "anti_periodic_condition",
            "non_attribuate_condition",
            "anti_periodic_condition",
            "homogeneous_Dirichlet_condition",
        ]
    ):
        # Initialize the boundary_condition list
        # Vertical boundary_condition
        boundary_condition_list[::N_point_theta] = 2
        boundary_condition_list[N_point_theta - 1 :: N_point_theta] = 3
        # Horizontal boundary_condition
        boundary_condition_list[N_point_theta * (N_point_r - 1) : N_total_point] = 1

        Periodic_point = -np.ones((N_point_r - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(
            0, N_point_theta * (N_point_r - 1), N_point_theta
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + N_point_theta - 1

    elif np.all(
        boundary_condition
        == [
            "anti_periodic_condition",
            "homogeneous_Dirichlet_condition",
            "anti_periodic_condition",
            "non_attribuate_condition",
        ]
    ):
        # Initialyze the boundary_condition list
        # Vertical boundary_condition
        boundary_condition_list[::N_point_theta] = 2
        boundary_condition_list[N_point_theta - 1 :: N_point_theta] = 3

        # Horizontal boundary_condition
        boundary_condition_list[:N_point_theta] = 1

        Periodic_point = -np.ones((N_point_r - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(
            N_point_theta, N_point_theta * (N_point_r), N_point_theta
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + N_point_theta - 1

    elif np.all(
        boundary_condition
        == [
            "periodic_condition",
            "non_attribuate_condition",
            "periodic_condition",
            "homogeneous_Dirichlet_condition",
        ]
    ):
        # Initialyze the boundary_condition list
        # Vertical boundary_condition
        boundary_condition_list[::N_point_theta] = 2
        boundary_condition_list[N_point_theta - 1 :: N_point_theta] = 2
        # Horizontal boundary_condition
        boundary_condition_list[N_point_theta * (N_point_r - 1) : N_total_point] = 1

        Periodic_point = -np.ones((N_point_r - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(
            0, N_point_theta * (N_point_r - 1), N_point_theta
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + N_point_theta - 1

    elif np.all(
        boundary_condition
        == [
            "periodic_condition",
            "homogeneous_Dirichlet_condition",
            "periodic_condition",
            "non_attribuate_condition",
        ]
    ):
        # Initialyze the boundary_condition list
        # Vertical boundary_condition
        boundary_condition_list[::N_point_theta] = 2
        boundary_condition_list[N_point_theta - 1 :: N_point_theta] = 2
        # Horizontal boundary_condition
        boundary_condition_list[:N_point_theta] = 1

        Periodic_point = -np.ones((N_point_r - 1, 2), dtype=np.int32)
        Periodic_point[:, 0] = np.arange(
            N_point_theta, N_point_theta * (N_point_r), N_point_theta
        )
        Periodic_point[:, 1] = Periodic_point[:, 0] + N_point_theta - 1

    else:
        raise NameError(
            "This boundary condition doesn't exist in this code or it's wrong"
        )

    return boundary_condition_list, Periodic_point
