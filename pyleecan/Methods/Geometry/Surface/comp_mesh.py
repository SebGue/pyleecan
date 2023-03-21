from math import ceil
from ....Functions.labels import decode_label, BOUNDARY_PROP_LAB, RADIUS_PROP_LAB


def comp_mesh(self, element_size, user_mesh_dict=None):
    """Returns the number of mesh elements on each line of the surface
    to match the element_size.

    Parameters
    ----------
    self : Surface
        a Surface object

    element_size : float
        The size of each element on the mesh [m]
    prefix: str
        Label prefix for line


    Returns
    -------
    mesh_list : list
        List containing the number of element of each line of the surface
    """
    mesh_dict = {}
    label = decode_label(self.label)["lam_type"]
    if user_mesh_dict and label in user_mesh_dict:
        mesh_dict.update(user_mesh_dict[label])

    mesh_list = []
    lines = self.get_lines()
    for line in lines:
        length = line.comp_length()
        number_of_element = ceil(length / element_size)
        if line.prop_dict and BOUNDARY_PROP_LAB in line.prop_dict:
            if line.prop_dict[BOUNDARY_PROP_LAB] in mesh_dict:
                number_of_element = mesh_dict[line.prop_dict[BOUNDARY_PROP_LAB]]
        elif line.prop_dict and RADIUS_PROP_LAB in line.prop_dict:
            if line.prop_dict[RADIUS_PROP_LAB] in mesh_dict:
                number_of_element = mesh_dict[line.prop_dict[RADIUS_PROP_LAB]]

        mesh_list.append(number_of_element)

    return mesh_list
