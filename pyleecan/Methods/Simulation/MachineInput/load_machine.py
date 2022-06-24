# -*- coding: utf-8 -*-

from pyleecan.Functions.load import load
import numpy as np


def load_machine(file_path):
    """Loads the geometrical input parameters of a PM machine

    Parameters
    ----------
    file_path: str
        path to the file of the machine to be loaded

    Returns
    -------
    active_length: float
        The active height of the PM machine
    height_magnet: float
        The height of the magnet of the PM machine
    width_magnet_av: float
        The width of the magnet of the PM machine
    stator_diameter_ext: float
        The outer diameter of the stator
    stator_diameter_int: float
        The inner diameter of the stator
    height_slotw: float
        The active height of the stator slot
    width_slotw_av: float
        The width of the stator slot, which is approximated as the ration of the slot surface to the slot height
    height_yoke_stator: float
        The height of the stator yoke
    stator_pole_pairs: float
        The number of the stator slots
    rotor_diameter_ext: float
        The outer diameter of the rotor
    rotor_diameter_int: float
        The inner diameter of the rotor
    rotor_pole_pairs:float
        The number of the rotor poles
    height_yoke_rotor: float
        The height of the rotor yoke
    """
    # Load Machine file
    Machine = load(file_path)
    """
    # Get the machine active length : supposed that active_length(rotor) = active_length(stator)
    active_length = Machine.rotor.L1

    # Get magnet geometrical parameters
    height_magnet = Machine.rotor.slot.comp_height_active()
    width_magnet_av = Machine.rotor.slot.comp_surface() / height_magnet

    # Get stator geometrical parameters
    stator_diameter_ext = Machine.stator.Rext
    stator_diameter_int = Machine.stator.Rint
    height_slotw = Machine.stator.slot.comp_height_active()
    width_slotw_av = Machine.stator.slot.comp_surface() / height_slotw

    height_yoke_stator = Machine.stator.comp_height_yoke()
    stator_pole_pairs = Machine.stator.get_Zs()

    # Get rotor geometrical parameters
    rotor_diameter_ext = Machine.rotor.Rext
    rotor_diameter_int = Machine.rotor.Rint
    rotor_pole_pairs = Machine.rotor.get_pole_pair_number()
    height_yoke_rotor = Machine.rotor.comp_height_yoke()

    e = Machine.stator.Rint - Machine.rotor.Rext

    return {
        "magnet_height": height_magnet,
        "magnet_average_width": width_magnet_av,
        "stator_diameter_ext": stator_diameter_ext,
        "stator_diameter_int": stator_diameter_int,
        "slot_height": height_slotw,
        "slot_width": width_slotw_av,
        "yoke_height_of_stator": height_yoke_stator,
        "slots_number": stator_pole_pairs,
        "rotor_diameter_ext": rotor_diameter_ext,
        "rotor_diameter_int": rotor_diameter_int,
        "rotor_pole_pairs": rotor_pole_pairs,
        "magnet_height": height_yoke_rotor,
        "e": e,
    }
    """
    tp = (
        2
        * np.tan(np.pi / Machine.rotor.get_pole_pair_number())
        * (Machine.rotor.Rint + 0.5 * Machine.rotor.comp_height_yoke())
    )

    # hm = 10e-3  # PM height in y direction (m)
    hm = Machine.rotor.slot.comp_height_active()

    # tm = 55e-3  # PM length in x direction (m)
    tm = (Machine.rotor.slot.comp_surface_active() / hm) * 0.000001

    # e = 1e-3  # Air-gap thickness (m)
    e = Machine.comp_width_airgap_mec()

    # hst = 30e-3  # Stator total height (m)
    hst = Machine.stator.comp_height_yoke()

    # hs = 20e-3  # Slot height (m)
    hs = Machine.stator.slot.comp_height_active()

    # hmbi = 10e-3  # Moving armature height (moving back iron height)
    hmbi = Machine.rotor.comp_height_yoke()

    # ws = 10e-3  # Slot opening (m)
    ws = Machine.stator.slot.comp_surface() / hs

    ts = 2 * ws  # Slot pitch (m)

    # la = 1  # Active length (m)
    la = Machine.rotor.L1

    return {tp, hm, tm, e, hst, hs, hmbi, ws, ts, la}