# -*- coding: utf-8 -*-
"""@package get_bore_line
@date Created on juin 20 10:56 2018
@author franco_i
"""
from numpy import pi, exp, linspace

from pyleecan.Classes.Arc1 import Arc1
from pyleecan.Classes.Arc3 import Arc3
from pyleecan.Methods import NotImplementedYetError


def get_bore_line(self, sym=1, label=""):
    """

    Parameters
    ----------
    self : Lamination
        a Lamination object
    alpha1 : float
        Starting angle [rad]
    alpha2 : float
        Ending angle [rad]
    label : str
        the label of the bore line

    Returns
    -------
    bore_line : list
        list of bore line

    """
    # adapted from LamSlotMulti.build_geometry
    bore_desc = self.get_bore_desc(sym=sym)
    bore_list = list()
    for bore in bore_desc:
        if type(bore["obj"]) is Arc1:
            bore_list.append(bore["obj"])
        else:
            lines = bore["obj"].build_geometry()
            for line in lines:
                line.rotate((bore["begin_angle"] + bore["end_angle"]) / 2)
            bore_list.extend(lines)

    return bore_list
