# -*- coding: utf-8 -*-

from numpy import sqrt


def comp_radius(self):
    """Compute the radius of the min and max circle that contains the slot

    Parameters
    ----------
    self : HoleM64
        A HoleM64 object

    Returns
    -------
    (Rmin,Rmax): tuple
        Radius of the circle that contains the slot [m]
    """
    Rbo = self.get_Rbo()
    Rmax = sqrt((Rbo - self.H0) ** 2 + self.H2**2)  # TODO include self.R0 (see HoleM60)
    Rmin = Rbo - (self.W1 + self.W2 + self.W3)

    return (Rmin, Rmax)
