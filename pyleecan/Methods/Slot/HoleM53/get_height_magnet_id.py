# -*- coding: utf-8 -*-


def get_height_magnet_id(self, index):
    """Get the height of the hole magnet of the corresponding index

    Parameters
    ----------
    self : HoleM53
        A HoleM53 object
    index : int
        Index of the magnet to get the height

    Returns
    -------
    Hmag: float
        Height of the Magnet [m]
    """

    # all magnets have the same height
    label = "magnet_" + str(index)
    if index in [0, 1] and getattr(self, label) is not None:
        return self.H2
    return 0
