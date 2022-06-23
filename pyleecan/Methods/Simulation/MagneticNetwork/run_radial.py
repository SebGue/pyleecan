
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:46:17 2022

@author: LAP02
"""
from solver_linear_model import linear_model
from geometry_linear_motor import geometry
from post_processing import compute_B_square,compute_B_radial
from plot import view_contour_flux,view_contour_flux2,Plot_B
import numpy as np
import meshio


la = 1            # Active length (m)
Br = 1.2          # PM remanent induction (residual induction) (T)
mu0 = np.pi*4e-7    # Permeability of vacuum (H/m)


pos=80
x_min = 0
x_max = 0.06


y_min = 0
y_max = 0.051

density=4
size_r= density*51+1
size_theta = density*60+1

r=np.linspace(0.005,0.05,size_r)
theta=np.linspace(0,np.pi/2,size_theta)

r_dual = (r[1:]+r[:-1])/2
theta_dual = (theta[1:]+theta[-1])/2

BC=["AP","HD","AP","HD"]
mode="polar"

F, list_geometry, Num_unknowns, list_elem, permeability_cell,list_coord = linear_model(
    size_theta, size_r,  theta,r, theta_dual,r_dual ,pos,BC,geometry,mu0,la,Br,mode)

x=(list_coord[:,1]*np.cos(list_coord[:,0])).reshape(size_r,size_theta)
y=(list_coord[:,1]*np.sin(list_coord[:,0])).reshape(size_r,size_theta)
list_cart=np.stack((x,y),axis=-1)

view_contour_flux(F,x,y,x.shape[1],x.shape[0], list_geometry)

Bx,By=compute_B_radial(F, list_elem,list_coord,la)

Norm=np.sqrt(Bx**2+By**2)


B=np.stack((Bx,By),axis=-1)

temp=np.zeros((list_elem.shape[0],3))
temp[:,0]=Bx.flatten()
temp[:,1]=By.flatten()
B=temp

temp=np.zeros((list_coord.shape[0],3))
temp[:,0:2]=np.array([x.flatten(),y.flatten()]).T

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
    "mymesh_radial.xdmf",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

# Alternative with the same options
#meshio.write_points_cells("mymesh.vtu", points, cells)

print("mesh saved",list_coord.shape,list_elem.shape)