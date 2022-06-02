# -*- coding: utf-8 -*-

from pyleecan.Functions.load import load_json


def load_parameters(file_path_json):
    """Load the stator and rotor parameters of an PM machine

    Parameters
    ----------
    file_path_json: str
        path to the file of the machine to be loaded

    Returns
    -------
    rotor_param: dict
        Geometrical parameters of the rotor core and magnets
    stator_param: dict
        Geometrical parameters of the stator core and stator slots
    """

    # Load the SPMSM machine from the .json file
    PM = load_json.load_json(file_path_json)

    # Get the dictionary values with the different motor parameters
    dict_PM = PM[1]

    # Get Rotor parameters
    rotor_param = {
        "rotor_Rext": dict_PM["rotor"]["Rext"],
        "rotor_Rint": dict_PM["rotor"]["Rint"],
        "rotor_poles_Zs": dict_PM["rotor"]["slot"]["Zs"],
        # Magnet height and width
        "magnet_height_Hmag": dict_PM["rotor"]["slot"]["Hmag"],
        "magnet_width_Wmag": dict_PM["rotor"]["slot"]["Wmag"],
    }

    # Get Stator parameters
    stator_param = {
        "stator_Rext": dict_PM["stator"]["Rext"],
        "stator_Rint": dict_PM["stator"]["Rint"],
        # Stator slot parameters
        "stator_tooth_tip_height_H0": dict_PM["stator"]["slot"]["H0"],
        "stator_tooth_tip_width_W0": dict_PM["stator"]["slot"]["W0"],
        "stator_slot_length_H2": dict_PM["stator"]["slot"]["H2"],
        "stator_slot_width_W2": dict_PM["stator"]["slot"]["W2"],
        "stator_poles_Zs": dict_PM["stator"]["slot"]["Zs"],
    }

    return rotor_param, stator_param
