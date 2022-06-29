from pyleecan.Classes.MagNetwork import MagNetwork
from pyleecan.Classes.MachineInput import load_machine
from pyleecan.Functions.load import load
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.InputCurrent import InputCurrent
from pyleecan.Classes.OPdq import OPdq
from pyleecan.Classes.ForceMT import ForceMT
from pyleecan.Functions.Plot import dict_2D
from pyleecan.Classes.MagFEMM import MagFEMM
from os.path import join
from Tests import save_validation_path as save_path


file_path = "C:/Users/pc/Downloads/SPMSM_val.json"

SPMSM_val = load(file_path)

load_machine(file_path)

# Create the simulation object of the MagNetwork
simu = Simu1(name="test_magnetwork", machine=SPMSM_val)

# Define Inputs and discretizations
# MagNetwork simulation
simu.input = InputCurrent(
    OP=OPdq(N0=1000, Id_ref=0, Iq_ref=0),
    Na_tot=60 * 4,
    Nt_tot=1,
)

# Definition of the magnetic simulation: with periodicity for the MagNetwork
simu.mag = MagNetwork(
    is_periodicity_a=True,
    is_periodicity_t=False,
)

# Mesh
# print(simu.mag.cartesianmeshclass_pyleecan())

# Test the MagNetwork simulation

simu.force = ForceMT()

# Create the simulation object of the MagNetwork
simu2 = Simu1(name="test_FEMM", machine=SPMSM_val)

# FEMM simulation
simu2.input = InputCurrent(
    OP=OPdq(N0=1000, Id_ref=0, Iq_ref=0),
    Na_tot=60 * 4,
    Nt_tot=1,
    is_periodicity_t=False,
    is_periodicity_a=True,
)
# Definition of the magnetic simulation: with periodicity for theFEMM
simu2.mag = MagFEMM(
    is_periodicity_a=True,
    is_periodicity_t=False,
    nb_worker=4,  # number of FEMM windows to be opened
    is_get_meshsolution=True,  # Get the mesh solution
    is_fast_draw=True,
    is_calc_torque_energy=False,
    Kmesh_fineness=2,  # Define the mesh fineness
)

simu2.force = ForceMT()

# Run both simulations
out = simu.run()  # MagNetwork
out2 = simu2.run()  # FEMM


# plot the 2D flux density curve
out.mag.B.plot_2D_Data("angle")
out2.mag.B.plot_2D_Data("angle")

# Compute the harmonics
"""
out.force.AGSF.plot_2D_Data(
    "wavenumber=[0,10]",
    "time[0]",
    data_list=[out.force.AGSF],
    legend_list=["Periodic", "Full"],
    is_show_fig=True,
    **dict_2D
)
"""

# Plot the results
out.mag.B.plot_2D_Data(
    # "time",
    "angle[0]{°}",
    # data_list=[out2.mag.B],
    data_list=[out.mag.B],
    legend_list=["Periodic", "Full"],
    is_show_fig=True,
    **dict_2D
)

out2.mag.B.plot_2D_Data(
    # "time",
    "angle[0]{°}",
    # data_list=[out2.mag.B],
    data_list=[out2.mag.B],
    legend_list=["Periodic", "Full"],
    save_path=join(save_path, simu2.name + "_B_space_FEMM.png"),
    is_show_fig=False,
    **dict_2D
)

# Compute the harmonics using the Maxwell tensor
out.force.AGSF.plot_2D_Data(
    "wavenumber=[0,10]",
    "time[0]",
    # data_list=[out2.force.AGSF],
    data_list=[out.force.AGSF],
    legend_list=["Periodic", "Full"],
    is_show_fig=True,
    **dict_2D
)

out2.force.AGSF.plot_2D_Data(
    "wavenumber=[0,10]",
    "time[0]",
    data_list=[out2.force.AGSF],
    legend_list=["Periodic", "Full"],
    save_path=join(save_path, simu.name + "_P_space_fft_FEMM.png"),
    is_show_fig=False,
    **dict_2D
)

# out.force.AGSF.plot_2D_Data(
#     "freqs",
#     "angle[0]",
#     data_list=[out2.force.AGSF],
#     legend_list=["Periodic", "Full"],
#     save_path=join(save_path, simu.name + "_B_space.png"),
#     is_show_fig=False,
#     **dict_2D
# )

pass
