# -*- coding: utf-8 -*-
from ....Classes.Segment import Segment
from ....Classes.Arc1 import Arc1


def build_geometry(self):
    """Compute the curve (Line) needed to plot the Slot.
    The ending point of a curve is the starting point of the next curve in
    the list

    Parameters
    ----------
    self : SlotM11
        A SlotM11 object

    Returns
    -------
    curve_list: list
        A list of 3 Segments

    """
    Rbo = self.get_Rbo()
    point_dict = self._comp_point_coordinate()
    Z1 = point_dict["Z1"]
    Z2 = point_dict["Z2"]
    Z3 = point_dict["Z3"]
    Z4 = point_dict["Z4"]
    ZM1 = point_dict["ZM1"]
    ZM2 = point_dict["ZM2"]
    ZM3 = point_dict["ZM3"]
    ZM4 = point_dict["ZM4"]

    # Creation of curve
    curve_list = []

    # Right side of slot
    if self.H0 > 0:
        if self.W1 == self.W0 and self.H1 < self.H0:
            curve_list.append(Segment(Z1, ZM2))
            curve_list.append(Segment(ZM2, Z2))
        else:
            curve_list.append(Segment(Z1, Z2))

    # Bottom of slot
    if self.is_outwards():
        if self.H0 > 0 and self.W1 != self.W0:
            curve_list.append(Arc1(Z2, ZM1, Rbo + self.H0, is_trigo_direction=True))
            curve_list.append(Arc1(ZM1, ZM4, Rbo + self.H0, is_trigo_direction=True))
            curve_list.append(Arc1(ZM4, Z3, Rbo + self.H0, is_trigo_direction=True))
        # curve_list.append(Arc1(Z2, Z3, Rbo - self.H0, is_trigo_direction=True))
        else:
            curve_list.append(Arc1(Z2, Z3, Rbo + self.H0, is_trigo_direction=True))
    else:
        if self.H0 > 0 and self.W1 != self.W0:
            curve_list.append(Arc1(Z2, ZM1, Rbo - self.H0, is_trigo_direction=True))
            curve_list.append(Arc1(ZM1, ZM4, Rbo - self.H0, is_trigo_direction=True))
            curve_list.append(Arc1(ZM4, Z3, Rbo - self.H0, is_trigo_direction=True))
        # curve_list.append(Arc1(Z2, Z3, Rbo - self.H0, is_trigo_direction=True))
        else:
            curve_list.append(Arc1(Z2, Z3, Rbo - self.H0, is_trigo_direction=True))

    # Left side of slot
    if self.H0 > 0:
        if self.W1 == self.W0 and self.H1 < self.H0:
            curve_list.append(Segment(Z3, ZM3))
            curve_list.append(Segment(ZM3, Z4))
        else:
            curve_list.append(Segment(Z3, Z4))

        # curve_list.append(Segment(Z3, Z4))

    return curve_list
