from ....Functions.Geometry.merge_notch_list import merge_notch_list


def get_magnet_list(self):
    """Returns an unordered list of magnets that are in the Lamination

    Parameters
    ----------
    self : Lamination
        A Lamination object

    Returns
    -------
    magnet_list : list
        list of magnets
    """
    magnet_list = []

    if hasattr(self, "magnet"):
        magnet_list.append(self.magnet)

    if hasattr(self, "hole"):
        for hole in self.hole:
            idx = 0
            while hasattr(hole, f"magnet_{idx}"):
                magnet = getattr(hole, f"magnet_{idx}")
                if magnet is not None:
                    magnet_list.append(magnet)
                idx += 1

    return magnet_list
