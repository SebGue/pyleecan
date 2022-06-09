# -*- coding: utf-8 -*-

from pyleecan.Functions.load import load


def load_machine(file_path):
    """Loads the stator and rotor material characteristics of a PM machine

    Parameters
    ----------
    file_path: str
        path to the file of the machine to be loaded

    Returns
    -------
    magnet_mat: float
        The remanence and permeability of the magnet
    rotor_permeance: ndarray/ float
        The non-linear (BH values) or linear material properties (mur_lin) of the rotor
    stator_permeance: ndarray/ float
        The non-linear (BH values) or linear material properties (mur_lin) of the stator
    winding_material: float/ str
        The permeability of the winding, if exists, else "Not defined" if it is not defind by the user
    insulation_material: float/ str
        The permeability of the insulation_material, if exists, else "Not defined" if it is not defind by the user
    shaft_mat: float/ str
        The non-linear (BH values) or linear material properties (mur_lin) of the shaft, if exists, else "Not defined" if it is not defind by the user
    """

    # Load the Machine from a defined path
    Machine = load(file_path)

    # Get the magnet material properties
    magnet_remanence = Machine.rotor.magnet.mat_type.mag.Brm20
    magnet_permeability = Machine.rotor.magnet.mat_type.mag.mur_lin

    # Get the rotor material properties
    if Machine.rotor.mat_type.is_isotropic == False:
        rotor_permeance = Machine.rotor.mat_type.mag.BH_curve.value
        # non-linear case
    else:
        rotor_permeance = Machine.rotor.mat_type.mag.mur_lin
        # linear case

    # Get the stator material properties
    if Machine.stator.mat_type.is_isotropic == False:
        stator_permeance = Machine.stator.mat_type.mag.BH_curve.value
    else:
        stator_permeance = Machine.stator.mat_type.mag.mur_lin

    # Get the conductor and insulator materials
    if Machine.stator.winding.conductor.cond_mat.mag != None:
        winding_material = Machine.stator.winding.conductor.cond_mat.mag.mur_lin
        insulation_material = Machine.stator.winding.conductor.ins_mat.mag.mur_lin
    else:
        winding_material = "Not defined"
        insulation_material = "Not defined"

    # Get the shaft material properties
    if Machine.shaft != None:
        if Machine.shaft.mat_type.is_isotropic == False:
            shaft_permeance = Machine.shaft.mat_type.mag.BH_curve.value
            # non-linear case
        else:
            shaft_permeance = Machine.shaft.mat_type.mag.mur_lin
            # linear case
    else:
        shaft_permeance = "Not defined"

    return {
        "magnet_remanence": magnet_remanence,
        "magnet_permeability": magnet_permeability,
        "rotor_permeance": rotor_permeance,
        "stator_permeance": stator_permeance,
        "winding_material": winding_material,
        "insulation_material": insulation_material,
        "shaft_permeance": shaft_permeance,
    }


print(load_machine("C:/Users/pc/Downloads/SPMSM_264s44p.json"))  # only for test