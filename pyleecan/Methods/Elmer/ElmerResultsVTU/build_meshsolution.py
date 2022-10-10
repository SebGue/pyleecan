# -*- coding: utf-8 -*-
from operator import is_
from numpy import arange, append as np_append
from meshio import read
from os.path import join, split, splitext

from SciDataTool import DataTime, Data1D, VectorField, Norm_ref

from ....Classes.SolutionData import SolutionData
from ....Classes.SolutionVector import SolutionVector
from ....Classes.MeshSolution import MeshSolution
from ....Classes.MeshVTK import MeshVTK

from ....Methods.Elmer.ElmerResultsVTU import ElmerResultsVTUError


# TODO add groups, see get_meshsolution of MagFEMM


def build_meshsolution(self, is_point_data=True):
    """Get the mesh and solution data from an Elmer VTU results file

    Parameters
    ----------
    self : ElmerResultsVTU
        a ElmerResultsVTU object

    is_point_data : bool
        True if the MeshSolution should be created from point data, else cell data

    Returns
    -------
    meshsol: MeshSolution
        MeshSolution created from the corresponding VTU file

    """
    # create meshsolution
    meshsol = MeshSolution(label=self.label)

    # get the mesh
    save_path, fn = split(self.file_path)
    file_name, file_ext = splitext(fn)
    if file_ext != ".vtu":
        raise ElmerResultsVTUError("ElmerResultsVTU: Results file must be of type VTU.")

    meshvtk = MeshVTK(path=save_path, name=file_name, format="vtu")
    # TODO maybe convert to MeshMat before
    meshsol.mesh = [meshvtk]

    # get the solution data on the mesh
    meshsolvtu = read(self.file_path)
    data = meshsolvtu.point_data if is_point_data else meshsolvtu.cell_data

    # setup axes
    if is_point_data:
        indices = arange(meshsolvtu.points.shape[0])
    else:
        indices = arange(
            meshsolvtu.cells[0].data.shape[0] + meshsolvtu.cells[1].data.shape[0]
        )
    Indices = Data1D(name="indice", values=indices, is_components=True)

    # store only data from store dict if available
    comp_ext = ["x", "y", "z"]

    sol_list = []  # list of solutions

    for key, value in data.items():
        # check if value should be stored
        if key in self.store_dict.keys():
            siz = value.shape[1] if is_point_data else value[0].shape[1]
            # only regard max. 3 components
            if siz > 3:
                self.get_logger().warning(
                    f'ElmerResultsVTU.build_meshsolution(): size of data "{key}" > 3'
                    + " - "
                    + "Data will be truncated."
                )
                siz = 3

            components = []
            comp_name = []
            if not is_point_data:
                value = np_append(value[0], value[1], axis=0)

            # loop though components
            for i in range(siz):
                # setup name, symbol and component name extension
                if siz == 1:
                    ext = ""
                else:
                    ext = comp_ext[i]

                # setup data object
                data = DataTime(
                    name=self.store_dict[key]["name"] + " " + ext,
                    unit=self.store_dict[key]["unit"],
                    symbol=self.store_dict[key]["symbol"] + ext,
                    axes=[Indices],
                    values=value[:, i],
                    normalizations={"ref": Norm_ref(ref=self.store_dict[key]["norm"])},
                )
                components.append(data)
                comp_name.append("comp_" + ext)

            # setup solution depending on number of field components
            if siz == 1:
                field = components[0]
                sol_list.append(
                    SolutionData(
                        field=field,
                        type_cell="point",
                        label=self.store_dict[key]["symbol"],
                    )
                )
            else:
                comps = {}
                for i in range(siz):
                    comps[comp_name[i]] = components[i]
                field = VectorField(
                    name=self.store_dict[key]["name"],
                    symbol=self.store_dict[key]["symbol"],
                    components=comps,
                )
                sol_list.append(
                    SolutionVector(
                        field=field,
                        type_cell="point",
                        label=self.store_dict[key]["symbol"],
                    )
                )

    meshsol.solution = sol_list

    return meshsol
