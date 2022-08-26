

import numpy as np
from pyleecan.Classes.MeshMat import MeshMat
from pyleecan.Classes.NodeMat import NodeMat
from pyleecan.Classes.CellMat import CellMat
from pyleecan.Classes.MeshSolution import MeshSolution
from pyleecan.Classes.SolutionMat import SolutionMat

def add_to_mesh_cell_math(self,Phi,Bx,By,list_elem,list_coord,type_coord_sys,material_dict):
    
    if type_coord_sys==2:
        x = list_coord[:, 1] * np.cos(list_coord[:, 0])
        y = list_coord[:, 1] * np.sin(list_coord[:, 0])
        list_coord[:, 0] = x
        list_coord[:, 1] = y

    mesh = MeshMat(dimension=3)
    mesh.node = NodeMat()
    print("Add points in NodeMath")
    for i in range(list_coord.shape[0]):
        mesh.node.add_node([list_coord[i,0],list_coord[i,1],0])
    print("Done \nAdd cells in MAthMesh")

    mesh.cell["quad"] = CellMat(nb_node_per_cell=4)
    for i in range(list_elem.shape[0]):
        mesh.add_cell(list_elem[i,:], "quad")

    print("Done \nAdd material for ech elementss")

    MSol = MeshSolution(mesh=[mesh])

    # for i in range(list_elem.shape[0]):
    #     MSol.group = {list_materials[list_geometry[i]-1]:list_elem[i,:]}

    print("Done")
    MSol.plot_mesh()

 #flux
    my_solution = SolutionMat(
        label="Flux (Weber)",
        type_cell="point",
        field=Phi,
        axis_name=[ "indice"],
        axis_size = [Phi.size],
    )
    MSol.solution.append(my_solution)
    MSol.plot_contour()

        #To do -> correct
    # B=np.zeros((list_elem.shape[0],2))
    # B[:,0]=Bx
    # B[:,1]=By

    # my_vec_solution = SolutionMat(
    #     label="my_vector",
    #     type_cell="cell", #type_cell="point"?
    #     field=B,
    #     incice=np.arange(list_elem.shape[0]),
    #     axis_name=["indice", "component"],
    #     axis_size = [B.shape],
    # )
    # MSol.solution.append(my_vec_solution)
    # MSol.plot_glyph(label="my_vector", is_point_arrow=True, factor=1/10)