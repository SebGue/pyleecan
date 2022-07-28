# from re import S
# import numpy as np
# import meshio
# import matplotlib.pyplot as plt
# from pyleecan.Functions.load import load

# from pyleecan.Classes.MeshMat import MeshMat
# from pyleecan.Classes.NodeMat import NodeMat
# from pyleecan.Classes.CellMat import CellMat
# from pyleecan.Classes.MeshSolution import MeshSolution
# from pyleecan.Classes.SolutionMat import SolutionMat

# def mesh_solution_pyleecan(self):

#     Launch mesh.io
#     list_coord = temp
#     points = list_coord
#     cells = [
#         ("quad", list_elem),
#     ]
#     mesh = meshio.Mesh(
#         points,
#         cells,
#         # Optionally provide extra data on points, cells, etc.
#         point_data={"Flux": Phi},
#         # Each item in cell data must match the cells array
#         cell_data={"B": [B], "Materials": [list_geometry]},
#     )
#     mesh.write(
#         "mymesh.xdmf",  # str, os.PathLike, or buffer/open file
#         # file_format="vtk",  # optional if first argument is a path; inferred from extension
#     )

#     Alternative with the same options
#     meshio.write_points_cells("mymesh.vtu", points, cells)


#     # Add my mesh to pyleecan
#     print("Solve RN done.")
#     mesh = MeshMat(dimension=3)
#     mesh.node = NodeMat()
#     print("Add points in mesh")
#     for i in range(list_cart.shape[0]):
#         mesh.node.add_node([list_cart[i, 0], list_cart[i, 1], 0])
#     print("Done \nAdd elements in mesh")

#     mesh.cell["quad"] = CellMat(nb_node_per_cell=4)
#     for i in range(list_elem.shape[0]):
#         mesh.add_cell(list_elem[i, :], "quad")

#     MSol = MeshSolution(mesh=[mesh])

#     # print("Done \nAdd material for ech elementss")
#     # for i in range(list_elem.shape[0]):
#     #     MSol.group = {list_materials[list_geometry[i] - 1]: list_elem[i, :]}

#     print("Done")
#     # plot the mesh
#     MSol.plot_mesh()

#     # plot the flux
#     field = Phi[np.newaxis]
#     print(field.shape)
#     my_solution = SolutionMat(
#         label="Flux (Weber)",
#         type_cell="point",
#         field=field,
#         indice=np.arange(list_coord.shape[0]),
#         axis_name=["time", "indice"],
#         axis_size=[1, list_coord.shape[0]],
#     )
#     MSol.solution.append(my_solution)
#     MSol.plot_contour()
