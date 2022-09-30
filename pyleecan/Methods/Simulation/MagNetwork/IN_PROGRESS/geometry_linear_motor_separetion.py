# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:51:29 2022

@author: LAP02
"""
from pre_processing import (
    save_mesh,
    init_reluc,
    init_point,
    init_cell,
    init_permeabilty_cell,
    init_permeabilty_cell,
    init_mesh_BC,
)
from assembler import assembly, right_member_assembly
import numpy as np
import matplotlib.pyplot as plt


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


#%%


list_materials = ["bob1", "bob2", "bob3", "air", "iron1", "PM", "iron3"]

permeability_materials = np.array([1, 1, 1, 1, 7500, 1, 7500])

x = 0.06
y = 0.051


def geometry(size_x, size_y, pos_pm):
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
    print(p0m)
    # Number of elements in the moving armature iron in y direction
    p0 = round(hmbi / h_y)
    # Number of elements in the magnetic air-gap (hm + e) in y direction
    p1 = round((hm + e) / h_y)

    m = p0 + p1
    n = size_x - 1
    nn = m * n
    print(nn)
    cells_materials = np.zeros(nn, dtype=np.uint16)

    mask_magnet = np.zeros(nn, dtype=np.bool_)
    mask_magnet[n * p1 - n : n * (p1 + p0 - p0m) + n] = True
    ### Geometry assembly
    for i in range(m):
        if p0 + p1 - p0m <= i < p0 + p1:
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

    reatach = np.arange(0, n + 1) + size_x * m
    return cells_materials, permeability_materials, reatach


def geometry2(size_x, size_y, pos_pm):
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
    print(p0m)
    # Number of elements in the moving armature iron in y direction
    p0 = round(hmbi / h_y)
    # Number of elements in the magnetic air-gap (hm + e) in y direction
    p1 = round((hm + e) / h_y)

    m = size_y - 1
    n = size_x - 1
    nn = (m - (p1 + p0)) * n
    print(nn)
    cells_materials = np.zeros(nn, dtype=np.uint16)

    mask_magnet = np.zeros(nn, dtype=np.bool_)
    mask_magnet[n * p1 - n : n * (p1 + p0 - p0m) + n] = True
    ### Geometry assembly
    i2 = 0
    for i in range(p1 + p0, m):

        if m - p0s <= i:
            for j in range(n):
                num_elem = n * i2 + j
                cells_materials[num_elem] = 5

        elif p0 + p1 <= i < n - p0s:
            for j in range(n):
                num_elem = n * i2 + j
                if m0s <= j < m0s + m1s:
                    cells_materials[num_elem] = 1
                elif 3 * m0s + m1s <= j < 3 * m0s + 2 * m1s:
                    cells_materials[num_elem] = 2
                elif 5 * m0s + 2 * m1s <= j < 5 * m0s + 3 * m1s:
                    cells_materials[num_elem] = 3
                else:
                    cells_materials[num_elem] = 5

        else:
            print("Wrong geometry")
        i2 += 1
    reatach = np.arange(0, n + 1)
    return cells_materials, permeability_materials, reatach


def numeroting_unknows(list_elem, BC_list, Periodic_point, reatach):
    """
    Pre-processing, numeroting the unknowns of the model in the linear system.

    Parameters
    ----------
    list_elem : nd-array, size: m (integers)
        tab of elements
    BC_list : nd-array, size: n (integers)
        data-strurture for BC evaluation.
    Periodic_point : nd-array, size: ? (int)
        Periodic point connection.

    Returns
    -------
    Num_Unknowns : nd-array, size: n (int)
        Numerotation in the linear system.

    """
    Num_Unknowns = -np.ones(BC_list.size, dtype=np.int32)
    mask = (BC_list != 1) & (BC_list != 3)
    mask[reatach] = False
    nn = mask.sum()
    Num_Unknowns[mask] = np.arange(nn, dtype=np.int32)

    Num_Unknowns[reatach] = np.arange(nn, nn + reatach.size, dtype=np.int32)

    Num_Unknowns[Periodic_point[:, 1]] = Num_Unknowns[Periodic_point[:, 0]]

    return Num_Unknowns


mode = "cartesian"
size_x = 61
size_y = 52
X = np.linspace(0, x, size_x)

######ROTOR###########

list_geometry_rotor, permeability_materials_rotor, reatach_rotor = geometry(
    size_x, size_y, 1
)
size_y_rotor = list_geometry_rotor.size // size_x + 2

plt.imshow(list_geometry_rotor.reshape((-1, size_x - 1)))
plt.show()

Y_rotor = np.linspace(0, hmbi + hm + e, size_y_rotor)
list_coord_rotor = init_point(size_x, size_y_rotor, X, Y_rotor)
list_elem_rotor = init_cell(size_x, size_y_rotor)
BC = ["AP", "HD", "AP", "NA"]
BC_list_rotor, Periodic_point_rotor = init_mesh_BC(size_x, size_y_rotor, BC)
Num_Unknowns_rotor = numeroting_unknows(
    list_elem_rotor, BC_list_rotor, Periodic_point_rotor, reatach_rotor
)

permeability_cell_rotor = init_permeabilty_cell(
    size_x, size_y_rotor, permeability_materials, mu0, list_geometry_rotor
)


reluc_list_rotor = init_reluc(list_elem_rotor, list_coord_rotor, mu0, la, mode)
M_csr_rotor = assembly(
    reluc_list_rotor,
    Num_Unknowns_rotor,
    list_elem_rotor,
    permeability_cell_rotor,
    BC_list_rotor,
)
print(list_coord_rotor[Periodic_point_rotor[:, 0]])
# print(list_coord_rotor[np.setdiff1d(np.arange(0,list_coord_rotor.shape[0]), reatach_rotor)])

#####STATOR##########

list_geometry_stator, permeability_materials_stator, reatach_stator = geometry2(
    size_x, size_y, 1
)
size_y_stator = list_geometry_stator.size // size_x + 2

plt.imshow(list_geometry_stator.reshape((-1, size_x - 1)))
plt.show()


Y_stator = np.linspace(hmbi + hm + e, y, size_y_stator)
list_coord_stator = init_point(size_x, size_y_stator, X, Y_stator)
list_elem_stator = init_cell(size_x, size_y_stator)
BC = ["AP", "NA", "AP", "HD"]
BC_list_stator, Periodic_point_stator = init_mesh_BC(size_x, size_y_stator, BC)

Num_Unknowns_stator = numeroting_unknows(
    list_elem_stator, BC_list_stator, Periodic_point_stator, reatach_stator
)
permeability_cell_stator = init_permeabilty_cell(
    size_x, size_y_stator, permeability_materials, mu0, list_geometry_stator
)


reluc_list_stator = init_reluc(list_elem_stator, list_coord_stator, mu0, la, mode)
M_csr_stator = assembly(
    reluc_list_stator,
    Num_Unknowns_stator,
    list_elem_stator,
    permeability_cell_stator,
    BC_list_stator,
)

print(list_coord_stator[Periodic_point_stator])

print(list_coord_stator[reatach_stator])
