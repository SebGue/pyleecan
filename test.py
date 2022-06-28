from pyleecan.Classes.MagNetwork import MagNetwork
from pyleecan.Classes.MachineInput import load_machine
from pyleecan.Functions.load import load
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.InputCurrent import InputCurrent
from pyleecan.Classes.OPdq import OPdq
from pyleecan.Classes.ForceMT import ForceMT
from pyleecan.Functions.Plot import dict_2D

file_path = "C:/Users/pc/Downloads/SPMSM_val.json"

SPMSM_val = load(file_path)

load_machine(file_path)

# Test the MagNetwork code
# Create simulation object
simu = Simu1(name="test_magnetwork", machine=SPMSM_val)

# Define Inputs and discretizations
simu.input = InputCurrent(
    OP=OPdq(N0=1000, Id_ref=0, Iq_ref=0),
    Na_tot=60 * 4,
    Nt_tot=1,
)

# Definition of the magnetic simulation: with periodicity
simu.mag = MagNetwork(
    is_periodicity_a=True,
    is_periodicity_t=False,
)
simu.force = ForceMT()

# Run simulation
out = simu.run()

# plot the 2D flux density curve
out.mag.B.plot_2D_Data("angle")

# Compute the harmonics
out.force.AGSF.plot_2D_Data(
    "wavenumber=[0,10]",
    "time[0]",
    data_list=[out.force.AGSF],
    legend_list=["Periodic", "Full"],
    is_show_fig=True,
    **dict_2D
)

pass
