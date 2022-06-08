# -*- coding: utf-8 -*-

from pyleecan.Functions.load import load


def load_machine(file_path):
    # Load Machine file
    Machine = load(file_path)

    # Get magnetic remanent flux density
    magnet_remanence = Machine.rotor.slot.material.get_Brm

    # Get stator material
    stator_material = Machine.stator.slot.material.get_BH

    # Get rotor material
    rotor_material = Machine.rotor.slot.material.get_BH

    return {
        "magnet_remanence": magnet_remanence,
        "stator_material": stator_material,
        "rotor_material": rotor_material,
    }


print(load_machine("C:/Users/pc/Downloads/SPMSM_264s44p.json"))  # only for test
