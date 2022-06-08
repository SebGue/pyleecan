# -*- coding: utf-8 -*-

from pyleecan.Functions.load import load


def load_machine(file_path):
    # Load Machine file
    Machine = load(file_path)

    # Get magnet geometrical parameters
    height_magnet = Machine.rotor.slot.comp_height_active()
    width_magnet_av = Machine.rotor.slot.comp_surface() / height_magnet

    # Get stator geometrical parameters
    height_slotw = Machine.stator.slot.comp_height_active()
    width_slotw_av = Machine.stator.slot.comp_surface() / height_slotw

    height_yoke_stator = Machine.stator.comp_height_yoke()
    stator_pole_pairs = Machine.stator.get_Zs()

    # Get rotor geometrical parameters
    rotor_pole_pairs = Machine.rotor.get_pole_pair_number()
    height_yoke_rotor = Machine.rotor.comp_height_yoke()

    return {
        "magnet_height": height_magnet,
        "magnet_average_width": width_magnet_av,
        "slot_height": height_slotw,
        "slot_width": width_slotw_av,
        "yoke_height_of_stator": height_yoke_stator,
        "slots_number": stator_pole_pairs,
        "rotor_pole_pairs": rotor_pole_pairs,
        "magnet_height": height_yoke_rotor,
    }


print(load_machine("C:/Users/pc/Downloads/SPMSM_264s44p.json"))  # only for test
