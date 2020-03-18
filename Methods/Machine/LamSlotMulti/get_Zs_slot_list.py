# -*- coding: utf-8 -*-

from numpy import pi, angle, exp

from pyleecan.Classes.Circle import Circle
from pyleecan.Classes.SurfLine import SurfLine
from pyleecan.Classes.Arc1 import Arc1
from pyleecan.Classes.Segment import Segment


def get_Zs_slot_list(self):
    """Return the number of slots of the individual slot in the list

    Parameters
    ----------
    self : LamSlotMulti
        a LamSlotMulti object

    Returns
    -------
    Zs_list : list
        Number of Slot

    """

    Zs_list = list()
    for slot in self.slot_list:
        Zs_list.append(slot.Zs)

    return Zs_list
