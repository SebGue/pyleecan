# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:46:17 2022

@author: LAP02
"""
from solver_linear_model import linear_model
from geometry_linear_motor import geometry
from post_processing import compute_B_square
import numpy as np
import meshio
from plot import view_contour_flux

la = 1            # Active length (m)
Br = 1.2          # PM remanent induction (residual induction) (T)
mu0 = np.pi*4e-7    # Permeability of vacuum (H/m)


pos=41
x_min = 0
x_max = 0.06


y_min = 0
y_max = 0.051

mul=2
size_x=int(60*mul)+1
size_y=int(51*mul)+1

x = np.linspace(x_min, x_max, size_x)
y = np.linspace(y_min, y_max, size_y)

x_dual = (x[1:]+x[:size_x-1])/2
y_dual = (y[1:]+y[:size_y-1])/2

BC=["P","HD","P","HD"]
mode="cartesian"

F, list_geometry, Num_unknowns, list_elem, permeability_cell,list_coord = linear_model(
    size_x, size_y, x, y, x_dual, y_dual, pos,BC,geometry,mu0,la,Br,mode)
print("Inconnues:",Num_unknowns.max()+1)
view_contour_flux(F,x,y,size_x,size_y, list_geometry)

Bx,By=compute_B_square(F, list_elem,list_coord,la)


temp=np.zeros((list_elem.shape[0],3))
temp[:,0]=Bx.flatten()
temp[:,1]=By.flatten()
B=temp

temp=np.zeros((list_coord.shape[0],3))
temp[:,0:2]=list_coord

list_coord=temp

points = list_coord
cells = [
    ("quad", list_elem),
]

mesh = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    point_data={"Flux": F},
    # Each item in cell data must match the cells array
    cell_data={"B": [B],"Materials":[list_geometry]},
)
mesh.write(
    "mymesh.xdmf",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

# Alternative with the same options
#meshio.write_points_cells("mymesh.vtu", points, cells)

print("mesh saved",list_coord.shape,list_elem.shape)