# -*- coding: utf-8 -*-

import numpy as np


def geometry_linear_motor(self):
    """
    Return the angular openings of the motor stator tooth, stator slot and rotor magnet.

    Parameters
    ----------
    self : MagNetwork
        A MagNetwork object
    """
    # Get machine object
    Machine = self.parent.machine

    # Compute the angular opening of the stator tooth
    angle_tooth = (
        2 * np.pi - Machine.stator.slot.comp_angle_opening() * Machine.stator.get_Zs()
    ) / Machine.stator.get_Zs()

    # Compute the angular opening of the stator slot
    angle_slot = Machine.stator.slot.comp_angle_opening()

    # Compute the angular opening of the rotor magnet
    angle_magnet = Machine.rotor.slot.comp_angle_opening()
    return {
        "angle_tooth": angle_tooth,
        "angle_slot": angle_slot,
        "angle_magnet": angle_magnet,
    }
