# -*- coding: utf-8 -*-

from pyleecan.Methods.Machine.Lamination.build_geometry import build_geometry
from pyleecan.Classes.MachineInput import load_machine
from pyleecan.Functions.load import load

file_path = "C:/Users/pc/Downloads/SPMSM_val.json"
SPMSM_val = load(file_path)

SPMSM_val.build_geometry(sym=2, alpha=0, delta=0, is_circular_radius=False)
