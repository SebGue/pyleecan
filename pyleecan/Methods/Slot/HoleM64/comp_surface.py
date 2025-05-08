# -*- coding: utf-8 -*-

from numpy import pi


def comp_surface(self):
    """Compute the Hole total surface (by analytical computation).

    Parameters
    ----------
    self : HoleM64
        A HoleM64 object

    Returns
    -------
    S: float
        Hole total surface [m**2]

    """
    # rectangle surface
    S1 = (self.W1 + self.W2 + self.W3) * self.H2
    # corner surface to substract
    S2 = (4 - pi) * self.R0**2

    return S1 - S2
