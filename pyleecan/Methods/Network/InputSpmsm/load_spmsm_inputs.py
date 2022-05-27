# -*- coding: utf-8 -*-

from pyleecan.Functions.load import load_json


def load_spmsm_inputs(file_path_json):
    """Load the stator and rotor characteristics of an SPMSM machine

    Parameters
    ----------
    file_path_json: str
        path to the file of the machine to be loaded

    Returns
    -------
    rotor_axial_length: float
        The axial length of the rotor [m]
    rotor_outer_radius: float
        The outer radius of the rotor [m]
    rotor_inner_radius: float
        The inner radius of the rotor [m]
    rotor_permeance_BH: float
        The B and H points of the BH material curve of the rotor (case of a non-linearity)
    rotor_permeance_linear: float
        The permeance value of the rotor (case of linearity)
    magnet_magnetization: int
        The magnetization direction of the magnet
    magnet_height: float
        The height of the magnet [m]
    magnet_width: float
        The width of the magnet [m]
    rotor_poles: int
        The number of the rotor poles
    rotor_poles: int
        The number of the rotor poles
    magnet_remanence: float
        The magnet remanence induction at 20degC (Default = 0) [T]
    magnet_permeability: float
        The relative magnetic permeability (Default = 1)
    stator_axial_length: float
        The axial length of the stator [m]
    stator_outer_radius: float
        The outer radius of the stator [m]
    stator_inner_radius: float
        The inner radius of the stator [m]
    stator_permeance_BH: float
        The B and H points of the BH material curve of the stator (case of a non-linearity)
    stator_permeance_linear: float
        The permeance value of the stator (case of linearity)
    stator_slots: int
        The number of the stator slots
    stator_slot_height: float
        The number of the stator slots
    stator_slot_angle: float
        The number of the stator slots
    """

    # Load the SPMSM machine
    SPMSM = load(file_path_json)

    return {
        # Rotor core properties for SPMSM (Status : done)
        "rotor_axial_length": SPMSM.rotor.L1,
        "rotor_outer_radius": SPMSM.rotor.Rext,
        "rotor_inner_radius": SPMSM.rotor.Rint,
        # Rotor core lamination material (Status : done)
        "rotor_permeance_BH": SPMSM.rotor.mat_type.mag.BH_curve.value,  # case non linear
        "rotor_permeance_linear": SPMSM.rotor.mat_type.mag.mur_lin,  # linear case
        # Magnets magnetization type (Status : done)
        "magnet_magnetization": SPMSM.rotor.magnet.type_magnetization,
        # Magnet height and width (Status : done)
        "magnet_height": SPMSM.rotor.slot.Hmag,
        "magnet_width": SPMSM.rotor.slot.Wmag,
        # Number of rotor poles (Status : done)
        "rotor_poles": SPMSM.rotor.slot.Zs,
        # Premanence and permeability of the magnet (Status : done)
        "magnet_remanence": SPMSM.rotor.magnet.mat_type.mag.Brm20,
        "magnet_permeability": SPMSM.rotor.magnet.mat_type.mag.mur_lin,
        # Stator core properties for SPMSM (Status : done)
        "stator_axial_length": SPMSM.stator.L1,
        "stator_outer_radius": SPMSM.stator.Rext,
        "stator_inner_radius": SPMSM.stator.Rint,
        # stator core lamination material (Status : done)
        "stator_permeance_BH": SPMSM.stator.mat_type.mag.BH_curve.value,  # case non linear
        "stator_permeance_linear": SPMSM.stator.mat_type.mag.mur_lin,  # linear case
        # Number of stator slots (Status : done)
        "stator_slots": SPMSM.stator.slot.Zs,
        # Parameters of the stator slots (Status : question to be asked :) )
        "stator_slot_height": SPMSM.stator.slot.H2,
        "stator_slot_angle": SPMSM.stator.slot.W2,
        # Parameters of the stator winding (Would they be useful?)
        "winding_layer": SPMSM.stator.winding.Nlayer,
        "parallel_circuits_per_phase": SPMSM.stator.winding.Npcp,
        "turns_per_coil": SPMSM.stator.Winding.Ntcoil,
        "coil_pitch": SPMSM.stator.Winding.coil_pitch,
        # Stator conductor (round wire)
        "stator_nb_wires": SPMSM.stator.Winding.conductor.Nwppc,
        "winding_diam_insulation": SPMSM.stator.Winding.conductor.Wins_cond,
        "winding_diam_no_insulation": SPMSM.stator.Winding.conductor.Wwire,
        "winding_strand_insulation": SPMSM.stator.Winding.conductor.Wins_wire,
        # Conductor and insulator materials
        "winding_material": SPMSM.stator.Winding.conductor.cond_mat.mag.mur_lin,
        "insulation_material": SPMSM.stator.Winding.conductor.ins_mat.mag.mur_lin,
    }
