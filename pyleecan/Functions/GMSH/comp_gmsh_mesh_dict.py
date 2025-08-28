from math import ceil, sqrt
from ...Functions.labels import BOUNDARY_PROP_LAB
from ...Functions.labels import short_label, decode_label


def comp_gmsh_mesh_dict(surface, element_size, user_mesh_dict={}):
    """Returns the number of mesh elements on each line of the surface
    to match the element_size.

    Parameters
    ----------
    self : Surface
        a Surface object

    element_size : float
        The default size of each element on the mesh [m]

    user_mesh_dict: dictionary
        User specified mesh properties


    Returns
    -------
    mesh_dict : dict
        Dictionary containing the number of element of each line of the surface
    """
    if user_mesh_dict is None:
        user_mesh_dict = {}

    mesh_dict = dict()
    lines = surface.get_lines()

    # get element size
    if surface.label in user_mesh_dict:
        element_size = user_mesh_dict[surface.label]

    # compute number of elements for all lines based on element size
    for ii, line in enumerate(lines):
        length = line.comp_length()
        number_of_element = ceil(length / element_size)
        mesh_dict[str(ii)] = number_of_element

        # overwrite number of elements given by boundary name in user_mesh_dict
        bc_prop = (
            line.prop_dict.get(BOUNDARY_PROP_LAB, None) if line.prop_dict else None
        )
        if bc_prop in user_mesh_dict:
            mesh_dict[str(ii)] = user_mesh_dict[bc_prop]

    return mesh_dict
