from numpy import exp, pi, cos, sin, tan, angle, sqrt


def _comp_point_coordinate(self):
    """Compute the point coordinates needed to plot the Slot.

    Parameters
    ----------
    self : HoleM65
        A HoleM65 object

    Returns
    -------
    point_dict: dict
        A dict of the slot coordinates
    """

    Rbo = self.get_Rbo()

    # comp point coordinate (in complex)

    # hole (on x axis)
    Z0 = -1j * (self.H1 / 2 - self.R)
    Z1 = self.R - 1j * (self.H1 / 2)
    Z2 = Z1.conjugate()
    Z3 = Z0.conjugate()
    Z4 = self.W2 - self.R - 1j * self.H1 / 2
    Z5 = self.W2 - 1j * (self.H1 / 2 - self.R)
    Z6 = Z5.conjugate()
    Z7 = Z4.conjugate()
    ZM1 = self.W2 - (self.W1 + self.W4) - 1j * self.H1 / 2
    ZM2 = self.W2 - self.W4 - 1j * self.H1 / 2
    ZM3 = ZM2.conjugate()
    ZM4 = ZM1.conjugate()

    # dZ the hole
    y0 = self.H1 / 2 + self.W3 / 2
    x0 = sqrt((Rbo - self.H0) ** 2 - y0**2)
    Z = x0 + 1j * y0
    dZ = x0 - self.W2 + 1j * y0

    Z0 += dZ
    Z1 += dZ
    Z2 += dZ
    Z3 += dZ
    Z4 += dZ
    Z5 += dZ
    Z6 += dZ
    Z7 += dZ
    ZM1 += dZ
    ZM2 += dZ
    ZM3 += dZ
    ZM4 += dZ

    # rotate to respect slot convention
    alpha = -pi / self.Zh
    Z *= exp(1j * alpha)
    Z0 *= exp(1j * alpha)
    Z1 *= exp(1j * alpha)
    Z2 *= exp(1j * alpha)
    Z3 *= exp(1j * alpha)
    Z4 *= exp(1j * alpha)
    Z5 *= exp(1j * alpha)
    Z6 *= exp(1j * alpha)
    Z7 *= exp(1j * alpha)
    ZM1 *= exp(1j * alpha)
    ZM2 *= exp(1j * alpha)
    ZM3 *= exp(1j * alpha)
    ZM4 *= exp(1j * alpha)

    # Draw the right hole by symmetry
    Zs = Z.conjugate()
    Z0s = Z0.conjugate()
    Z1s = Z1.conjugate()
    Z2s = Z2.conjugate()
    Z3s = Z3.conjugate()
    Z4s = Z4.conjugate()
    Z5s = Z5.conjugate()
    Z6s = Z6.conjugate()
    Z7s = Z7.conjugate()
    ZM1s = ZM1.conjugate()
    ZM2s = ZM2.conjugate()
    ZM3s = ZM3.conjugate()
    ZM4s = ZM4.conjugate()

    point_dict = dict()
    point_dict["Z"] = Z
    point_dict["Z0"] = Z0
    point_dict["Z1"] = Z1
    point_dict["Z2"] = Z2
    point_dict["Z3"] = Z3
    point_dict["Z4"] = Z4
    point_dict["Z5"] = Z5
    point_dict["Z6"] = Z6
    point_dict["Z7"] = Z7
    point_dict["ZM1"] = ZM1
    point_dict["ZM2"] = ZM2
    point_dict["ZM3"] = ZM3
    point_dict["ZM4"] = ZM4

    point_dict["Zs"] = Zs
    point_dict["Z0s"] = Z0s
    point_dict["Z1s"] = Z1s
    point_dict["Z2s"] = Z2s
    point_dict["Z3s"] = Z3s
    point_dict["Z4s"] = Z4s
    point_dict["Z5s"] = Z5s
    point_dict["Z6s"] = Z6s
    point_dict["Z7s"] = Z7s
    point_dict["ZM1s"] = ZM1s
    point_dict["ZM2s"] = ZM2s
    point_dict["ZM3s"] = ZM3s
    point_dict["ZM4s"] = ZM4s

    return point_dict
