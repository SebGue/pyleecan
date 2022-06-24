from pyleecan.Classes.MagneticNetwork import MagneticNetwork

from pyleecan.Methods.Simulation.MagneticNetwork import solver_linear_model
from pyleecan.Methods.Simulation.MagneticNetwork import run_radial

file_path = "C:/Users/pc/Downloads/SPMSM_val.json"
result = MagneticNetwork()
result.run_radial()
