# -*- coding: utf-8 -*-

from pyleecan.Methods.Machine.Lamination.build_geometry import build_geometry
from pyleecan.Functions.load import load
from pyleecan.Classes.LamSlotWind import LamSlotWind

file_path = "C:/Users/pc/Downloads/SPMSM_val.json"
machine = load(file_path)

# Machine get_lam_list method
machine.get_lam_list()

# Sorted stator and rotor lists
machine.get_lam_list(key="Stator")
machine.get_lam_list(key="Rotor")

# Get machine labels
machine.get_lam_list_label()


# machine.build_geometry(sym=2, alpha=0, delta=0, is_circular_radius=False)

machine.build_geometry(sym=2, alpha=0, delta=0)