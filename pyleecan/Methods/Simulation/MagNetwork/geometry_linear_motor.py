# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:51:29 2022

@author: LAP02
"""

import numpy as np
from pyleecan.Functions.load import load


def geometry_linear_motor(self, size_x, size_y, pos_pm):

    # Get machine object
    Machine = self.parent.machine

    # Definition of the machine geometrical input parameters
    # tp = 60e-3  # pole pitch (m)
    tp = (
        np.pi * Machine.stator.Rint / Machine.rotor.get_pole_pair_number()
    )  # Ref:https://www.slideshare.net/monprado1/11-basic-concepts-of-a-machine-77442134

    # tp = 2 * np.pi / Machine.comp_periodicity_spatial()[0]

    # hm = 10e-3  # PM height in y direction (m)
    hm = Machine.rotor.slot.comp_height_active()

    # tm = 55e-3  # PM width in x direction (m)
    tm = 2 * Machine.rotor.slot.comp_surface_active() / hm
    #  tm = Machine.rotor.slot.comp_width_opening(is_curved=True)

    # e = 1e-3  # Air-gap thickness (m)
    e = Machine.comp_width_airgap_mec()

    # hs = 20e-3  # Stator slot height (m)
    hs = Machine.stator.slot.comp_height()

    # hst = 30e-3  # Stator total height (m)
    hst = Machine.stator.Rext - Machine.stator.Rint

    # hmbi = 10e-3  # Moving armature height (moving back iron height)
    hmbi = Machine.rotor.comp_height_yoke()

    # ws = 10e-3  # Slot opening (m)
    ws = Machine.stator.slot.comp_width()

    # Stator yoke
    sy = Machine.stator.comp_height_yoke()

    # Stator tooth width (m)
    wt = (
        Machine.stator.get_Rbo()
        * 2
        * np.sin(
            2 * np.pi / Machine.stator.get_Zs()
            - float(np.arcsin(ws / (2 * Machine.stator.get_Rbo())))
        )
    )

    ts = ws + wt  # Slot pitch (m)

    # la = 1  # Active length (m)
    la = Machine.rotor.L1

    # Material properties
    # Br = 1.2  # PM remanent induction (residual induction) (T)
    Br = Machine.rotor.magnet.mat_type.mag.Brm20

    mu0 = np.pi * 4e-7  # Permeability of vacuum (H/m)
    mur1 = 1  # Relative permeability of air

    # mur2 = 7500 # Relative permeability of the stator iron
    mur2 = Machine.stator.mat_type.mag.mur_lin

    # Relative permeability of the rotor iron
    # mur3 = 7500
    mur3 = Machine.rotor.mat_type.mag.mur_lin

    # Relative permeability of the winding
    if Machine.stator.winding.conductor.cond_mat.mag != None:
        mur_bob = Machine.stator.winding.conductor.cond_mat.mag.mur_lin
    else:
        mur_bob = 1

    # Relative permeabiltity of the PM
    mur_PM = Machine.rotor.magnet.mat_type.mag.mur_lin

    # Generalize the materials list and their linear permeabilities for a specific number of windings
    list_materials = []
    permeabilty_materials = np.array([])

    for i in range(
        (int)(Machine.stator.get_Zs() / (2 * Machine.comp_periodicity_spatial()[0]))
    ):
        # list_materials = ["bob1", "bob2", "bob3", "air", "iron1", "PM", "iron3"]
        list_materials.append("bob[i]")
        np.append(permeabilty_materials, [mur_bob])

    list_materials += ["air", "iron1", "PM", "iron3"]

    np.append(permeabilty_materials, [mur1, mur2, mur_PM, mur3])

    # x and y positions
    x = tp
    y = Machine.stator.Rext - Machine.rotor.Rint

    # Definition of x-axis and y-axis steps
    # h_x = x / (size_x - 1)
    h_theta = (2 * np.pi / Machine.comp_periodicity_spatial()[0]) / (size_x - 1)
    h_y = y / (size_y - 1)

    # Number of elements in the stator armature
    # Number of elements in 1/2 tooth in x-axis direction
    # m0s = round(wt / 2 / h_x)
    # m1s = round(wt / h_x)
    angle_t = self.generalize_geometry[
        "angle_tooth"
    ]  # output from the generalize_geometry method to calculate the stator tooth opening
    theta_m0s = round(angle_t / 2 / h_theta)
    theta_m1s = round(angle_t / h_theta)

    # Number of elements in the stator back-iron in y-axis direction, generalized
    p0s = round(sy / h_y)

    # Total number of elements of the stator in x direction
    m = Machine.stator.get_Zs() * theta_m0s

    # Total number of element of stator in y direction
    p = (int)(
        Machine.stator.get_Zs() / (2 * Machine.rotor.comp_periodicity_spatial()[0])
    ) * p0s

    # Number of elements in the moving armature (the air-gap is supposed to be part of the moving armature)
    # Number of elements in half the air-gap between two adjacent PM in x direction
    m0m = round((tp - tm) / 2 / h_theta)
    m1m = round(tm / 2 / h_theta)
    p0m = round(e / h_y)  # Number of elements in the air-gap in y direction

    # Number of elements in the moving armature iron in y direction
    p0 = round(hmbi / h_y)

    # Number of elements in the magnetic air-gap (hm + e) in y direction
    p1 = round((hm + e) / h_y)

    # Number of elements in the y-axis direction
    m = size_y - 1

    # Number of elements in the x-axis direction
    n = size_x - 1

    # Total number of elements
    nn = m * n

    # Attribute an array of nn size of zeros to cells_materials
    cells_materials = np.zeros(nn, dtype=np.uint16)

    # Attribute an array of nn size of False statements to mask_magnet
    mask_magnet = np.zeros(nn, dtype=np.bool_)

    # Assign true entities tothe elements representing the
    mask_magnet[n * p1 - n : n * (p1 + p0 - p0m) + n] = True

    ### Geometry assembly
    for i in range(m):
        if m - p0s <= i:
            for j in range(n):
                num_elem = n * i + j
                cells_materials[num_elem] = 5

        elif p0 + p1 <= i < n - p0s:
            for j in range(n):
                num_elem = n * i + j
                if theta_m0s <= j < theta_m0s + theta_m1s:
                    cells_materials[num_elem] = 1
                elif 3 * theta_m0s + theta_m1s <= j < 3 * theta_m0s + 2 * theta_m1s:
                    cells_materials[num_elem] = 2
                elif 5 * theta_m0s + 2 * theta_m1s <= j < 5 * theta_m0s + 3 * theta_m1s:
                    cells_materials[num_elem] = 3
                else:
                    cells_materials[num_elem] = 5

        elif p0 + p1 - p0m <= i < p0 + p1:
            for j in range(n):
                num_elem = n * i + j
                cells_materials[num_elem] = 4
            ##
        elif p1 <= i < p1 + p0 - p0m:
            for j in range(n):
                num_elem = n * i + j
                if pos_pm + 2 * m1m >= n:
                    if (pos_pm + 2 * m1m) % n < j <= pos_pm % n:
                        cells_materials[num_elem] = 4
                    else:
                        cells_materials[num_elem] = 6
                else:
                    if pos_pm <= j < (pos_pm + 2 * m1m):
                        cells_materials[num_elem] = 6
                    else:
                        cells_materials[num_elem] = 4
            ##
        elif i < p1:
            for j in range(n):
                num_elem = n * i + j
                cells_materials[num_elem] = 7
        else:
            print("Wrong geometry")

    return cells_materials, permeabilty_materials
