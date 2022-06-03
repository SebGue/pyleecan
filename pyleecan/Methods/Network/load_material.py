# -*- coding: utf-8 -*-

from pyleecan.Functions.load import load_json


def LoadMaterial(file_path_json):
    """Load the stator and rotor material characteristics of a PM machine

    Parameters
    ----------
    file_path_json: str
        path to the file of the machine to be loaded

    Returns
    -------
    magnet_mat: dict
        The magnet_remanence and magnet_permeability of the magnet
    rotor_mat: dict
        The linear or non-linear material properties of the rotor
    stator_mat: dict
        The linear or non-linear material properties of the stator
    winding_mat: dict
        The material properties of the winding and insulation
    shaft_mat: dict
        The linear or non-linear material properties of the shaft
    """

    # Load the PM machine from the .json file
    PM = load_json.load_json(file_path_json)

    # Get the dictionary values with the different motor parameters
    dict_PM = PM[1]

    # Get Magnet material properties
    magnet_mat = {
        "magnet_remanence": dict_PM["rotor"]["magnet"]["mat_type"]["mag"]["Brm20"],
        "magnet_permeability": dict_PM["rotor"]["magnet"]["mat_type"]["mag"]["mur_lin"],
    }

    # Get Rotor material properties
    if dict_PM["rotor"]["mat_type"]["is_isotropic"] == False:
        rotor_mat = {
            "rotor_permeance": dict_PM["rotor"]["mat_type"]["mag"]["BH_curve"]["value"]
        }  # non-linear case
    else:
        rotor_mat = {
            "rotor_permeance": dict_PM["rotor"]["mat_type"]["mag"]["mur_lin"]
        }  # linear case

    # Get Stator material properties
    if dict_PM["stator"]["mat_type"]["is_isotropic"] == False:
        stator_mat = {
            "stator_permeance": dict_PM["stator"]["mat_type"]["mag"]["BH_curve"][
                "value"
            ]
        }
    else:
        stator_mat = {
            "stator_permeance": dict_PM["stator"]["mat_type"]["mag"]["mur_lin"]
        }

    # Conductor and insulator materials
    if dict_PM["stator"]["winding"]["conductor"]["cond_mat"]["mag"] != None:
        winding_mat = {
            "winding_material": dict_PM["stator"]["winding"]["conductor"]["cond_mat"][
                "mag"
            ]["mur_lin"],
            "insulation_material": dict_PM["stator"]["winding"]["conductor"]["ins_mat"][
                "mag"
            ]["mur_lin"],
        }
    else:
        winding_mat = {"winding_material": "Not defined"}

    # Get Shaft material properties
    if dict_PM["shaft"] != None:
        if dict_PM["shaft"]["mat_type"]["is_isotropic"] == False:
            shaft_mat = {
                "shaft_permeance": dict_PM["shaft"]["mat_type"]["mag"]["BH_curve"][
                    "value"
                ]
            }  # non-linear case
        else:
            shaft_mat = {
                "shaft_permeance": dict_PM["shaft"]["mat_type"]["mag"]["mur_lin"]
            }  # linear case
    else:
        shaft_mat = {"shaft_permeance": "Not defined"}

    return magnet_mat, rotor_mat, stator_mat, winding_mat, shaft_mat