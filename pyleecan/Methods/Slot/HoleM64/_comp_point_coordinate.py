from numpy import exp, pi, cos, sin, tan, angle, sqrt


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

    # comp point coordinate (in complex)
    # 4 corner points of the complete hole
    Zc1 = self.Rbo - self.H0 - 1j * self.W2
    Zc2 = self.Rbo - self.H0 + 1j * self.W2
    Zc3 = self.Rbo - self.H0 - W + 1j * self.W2
    Zc4 = self.Rbo - self.H0 - W - 1j * self.W2

    # include the radii
    Z11 = Zc1 - self.R0
    Z12 = Zc1 + 1j * self.R0

    Z21 = Zc2 - 1j * self.R0
    Z22 = Zc2 - self.R0

    Z31 = Zc3 + self.R0
    Z32 = Zc3 - 1j * self.R0

    Z41 = Zc4 + 1j * self.R0
    Z42 = Zc4 + self.R0

    # 4 corner points of the magnet
    ZM1 = Zc1 - self.W1
    ZM2 = Zc2 - self.W1
    ZM3 = ZM2 - self.W2
    ZM4 = ZM1 - self.W2

    # Rotation of angle -pi / 2 to respect slot conventions
    rotation = exp(-1j * pi / 2)
    Z11 *= rotation
    Z12 *= rotation
    Z21 *= rotation
    Z22 *= rotation
    Z31 *= rotation
    Z32 *= rotation
    Z41 *= rotation
    ZM42 *= rotation
    ZM1 *= rotation
    ZM2 *= rotation
    ZM3 *= rotation
    ZM4 *= rotation

    point_dict = dict()
    point_dict["Z11"] = Z11
    point_dict["Z12"] = Z12
    point_dict["Z21"] = Z21
    point_dict["Z22"] = Z22
    point_dict["Z31"] = Z31
    point_dict["Z32"] = Z32
    point_dict["Z41"] = Z41
    point_dict["Z42"] = Z42
    point_dict["ZM1"] = ZM1
    point_dict["ZM2"] = ZM2
    point_dict["ZM3"] = ZM3
    point_dict["ZM4"] = ZM4

    return point_dict
