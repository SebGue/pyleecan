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
    # 4 corner points of the complete hole
    Zc1 = x0
    Zc2 = x0 + 1j * self.H2 / 2
    Zc3 = x0 - W + 1j * self.H2 / 2
    Zc4 = x0 - W

    # include the radii
    Z1 = Zc1

    Z21 = Zc2 - 1j * self.R0
    Z22 = Zc2 - self.R0
    Z2c = Z21 - self.R0

    Z31 = Zc3 + self.R0
    Z32 = Zc3 - 1j * self.R0
    Z3c = Z32 + self.R0

    Z4 = Zc4

    # 4 corner points of the magnet
    ZM1 = Zc1 - self.W1
    ZM2 = Zc2 - self.W1
    ZM3 = ZM2 - self.W2
    ZM4 = ZM1 - self.W2

    point_dict = dict()
    point_dict["Z1"] = Z1
    point_dict["Z21"] = Z21
    point_dict["Z22"] = Z22
    point_dict["Z2c"] = Z2c
    point_dict["Z31"] = Z31
    point_dict["Z32"] = Z32
    point_dict["Z3c"] = Z3c
    point_dict["Z4"] = Z4
    point_dict["ZM1"] = ZM1
    point_dict["ZM2"] = ZM2
    point_dict["ZM3"] = ZM3
    point_dict["ZM4"] = ZM4

    rotate = 1  # exp(-1j*pi/self.Zh)
    for key in point_dict.keys():
        point_dict[key] *= rotate

    return point_dict
