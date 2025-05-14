from numpy import pi, exp


def _comp_point_coordinate(self):
    """Compute the point coordinates needed to plot the Slot.

    Parameters
    ----------
    self : HoleM64
        A HoleM64 object

    Returns
    -------
    point_dict: dict
        A dict of the slot coordinates
    """

    Rbo = self.get_Rbo()
    W = self.W1 + self.W2 + self.W3
    x0 = Rbo - self.H0

    # comp point coordinate (in complex)
    # 2 corner points of the complete hole
    Zc1 = x0 + 1j * self.H2 / 2
    Zc2 = x0 - W + 1j * self.H2 / 2

    # include the radii
    Z11 = Zc1 - 1j * self.R0
    Z12 = Zc1 - self.R0
    Z1c = Z11 - self.R0

    Z21 = Zc2 + self.R0
    Z22 = Zc2 - 1j * self.R0
    Z2c = Z22 + self.R0

    # 4 corner points of the magnet
    ZM1 = Zc1 - self.W1
    ZM2 = ZM1 - self.W2

    point_dict = dict()
    point_dict["Z11"] = Z11
    point_dict["Z12"] = Z12
    point_dict["Z1c"] = Z1c

    point_dict["Z21"] = Z21
    point_dict["Z22"] = Z22
    point_dict["Z2c"] = Z2c

    point_dict["ZM1"] = ZM1
    point_dict["ZM2"] = ZM2

    rotate = 1  # exp(-1j*pi/self.Zh)
    for key in point_dict.keys():
        point_dict[key] *= rotate

    return point_dict
