from pyleecan.Classes.MagNetwork import MagNetwork
from pyleecan.Classes.MachineInput import load_machine

file_path = "C:/Users/pc/Downloads/SPMSM_val.json"

load_machine(file_path)

MagNetwork(file_path)

MagNetwork(file_path).geometry_linear_motor
# MagNetwork(file_path).solver_linear_model

MagNetwork(file_path).run_radial()