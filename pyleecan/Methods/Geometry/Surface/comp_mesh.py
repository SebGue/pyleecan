from math import ceil


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
    mesh_dict : dict
        Dictionary containing the number of element of each line of the surface
    """
    mesh_dict = {}
    if user_mesh_dict and self.label in user_mesh_dict:
        mesh_dict.update(user_mesh_dict[self.label])

    mesh_list = []
    lines = self.get_lines()
    for line in lines:
        length = line.comp_length()
        number_of_element = ceil(length / element_size)
        mesh_list.append(number_of_element)

    return mesh_list
