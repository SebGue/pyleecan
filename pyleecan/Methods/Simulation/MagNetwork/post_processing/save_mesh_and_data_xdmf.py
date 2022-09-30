<<<<<<< HEAD
=======

>>>>>>> c00deaf7caddf7ae6eba9d88d694c0b0f47f0a8a
import numpy as np
import meshio


<<<<<<< HEAD
def save_mesh_and_data_xdmf(
    self, file_name, Phi, Bx, By, list_elem, list_coord, type_coord_sys, material_dict
):
    print("Launch save xdmf")
    if type_coord_sys == 2:
=======

def save_mesh_and_data_xdmf(self,file_name,Phi,Bx,By,list_elem,list_coord,type_coord_sys,material_dict):
    print("Launch save xdmf")
    if type_coord_sys==2:
>>>>>>> c00deaf7caddf7ae6eba9d88d694c0b0f47f0a8a
        x = list_coord[:, 1] * np.cos(list_coord[:, 0])
        y = list_coord[:, 1] * np.sin(list_coord[:, 0])
        list_coord[:, 0] = x
        list_coord[:, 1] = y

<<<<<<< HEAD
    B = np.zeros((list_elem.shape[0], 3))
    B[:, 0] = Bx
    B[:, 1] = By

    points = list_coord
    # Define Mesh.io cells
=======
    B=np.zeros((list_elem.shape[0],3))
    B[:,0]=Bx
    B[:,1]=By

    
    points = list_coord
    #Define Mesh.io cells
>>>>>>> c00deaf7caddf7ae6eba9d88d694c0b0f47f0a8a
    cells = [
        ("quad", list_elem),
    ]

<<<<<<< HEAD
    mag = np.sqrt(Bx ** 2 + By ** 2)
    mesh = meshio.Mesh(
        points,
        cells,
        # Optionally provide extra data on points, cells, etc.
        point_data={"Flux": Phi},
        # Each item in cell data must match the cells array
        # cell_data={"B": [B],"Materials":[materials_dict],"Permability":[list_elem_permability]},
        cell_data={"B": [B], "B_magnitude": [mag]},
    )
    mesh.write(
        file_name + ".xdmf",
    )
    print("Done, The xdmf is avaible.")
    return
=======
    mag=np.sqrt(Bx**2+By**2)
    mesh = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    point_data={"Flux": Phi},
    # Each item in cell data must match the cells array
    #cell_data={"B": [B],"Materials":[materials_dict],"Permability":[list_elem_permability]},
    
    cell_data={"B": [B],"B_magnitude":[mag]}
    )
    mesh.write(
        file_name+".xdmf",  
    )
    print("Done, The xdmf is avaible.")
    return 
>>>>>>> c00deaf7caddf7ae6eba9d88d694c0b0f47f0a8a
