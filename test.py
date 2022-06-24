from pyleecan.Classes.MagneticNetwork import MagneticNetwork
from pyleecan.Classes.MachineInput import load_machine

file_path = "C:/Users/pc/Downloads/SPMSM_val.json"

load_machine(file_path)

MagneticNetwork(file_path)

MagneticNetwork(file_path).geometry_linear_motor

MagneticNetwork(file_path).run_radial