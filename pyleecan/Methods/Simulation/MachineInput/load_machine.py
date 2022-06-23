# -*- coding: utf-8 -*-

from pyleecan.Functions.load import load


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
    }
