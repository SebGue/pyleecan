# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:51:29 2022
@author: LAP02
"""

import numpy as np

import numpy as np


def geometry_motor(self, size_x, size_y, pos_pm):

    # Machine's characteristics
    tp = 60e-3  # pole pitch (m)
    tm = 55e-3  # PM length in x direction (m)
    hm = 10e-3  # PM height in y direction (m)
    e = 1e-3  # Air-gap thickness (m)
    hst = 30e-3  # Stator total height (m)
    hs = 20e-3  # Slot height (m)
    hmbi = 10e-3  # Moving armature height (moving back iron height)
    ws = 10e-3  # Slot opening (m)
    ts = 2 * ws  # Slot pitch (m)
    la = 1  # Active length (m)
    Br = 1.2  # PM remanent induction (residual induction) (T)
    mu0 = np.pi * 4e-7  # Permeability of vacuum (H/m)
    mur1 = 1  # Relative permeability of air
    mur2 = 7500  # Relative permeability of statot iron
    mur3 = 7500

    list_materials = ["bob1", "bob2", "bob3", "air", "iron1", "PM", "iron3"]

    permeabilty_materials = np.array([1, 1, 1, 1, 7500, 1, 7500])

    x = 0.06
    y = 0.051

    print(size_x, size_y)
    h_x = x / (size_x - 1)
    h_y = y / (size_y - 1)
    # % Number of elements in the stator armature
    m0s = round(
        (ts - ws) / 2 / h_x
    )  # Number of elements in half a tooth in x direction
    m1s = round(ws / h_x)
    # Number of elements in the stator back iron in y direction
    p0s = round((hst - hs) / h_y)
    m = 12 * m0s  # Total number of elements of the stator in x direction
    p = 3 * p0s  # Total number of element of stator in y direction

    # Number of elements in the moving armature (the air-gap is supposed to be part of the moving armature)
    # Number of elements in half the air-gap between two adjacent PM in x direction
    m0m = round((tp - tm) / 2 / h_x)
    m1m = round(tm / 2 / h_x)
    p0m = round(e / h_y)  # Number of elements in the air-gap in y direction
    # Number of elements in the moving armature iron in y direction
    p0 = round(hmbi / h_y)
    # Number of elements in the magnetic air-gap (hm + e) in y direction
    p1 = round((hm + e) / h_y)

    m = size_y - 1
    n = size_x - 1
    nn = m * n
    print(nn)
    cells_materials = np.zeros(nn, dtype=np.uint16)

    mask_magnet = np.zeros(nn, dtype=np.bool_)
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
                if m0s <= j < m0s + m1s:
                    cells_materials[num_elem] = 1
                elif 3 * m0s + m1s <= j < 3 * m0s + 2 * m1s:
                    cells_materials[num_elem] = 2
                elif 5 * m0s + 2 * m1s <= j < 5 * m0s + 3 * m1s:
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
