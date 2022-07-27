from pyleecan.Classes.MagNetwork import MagNetwork
from pyleecan.Functions.load import load
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.InputCurrent import InputCurrent
from pyleecan.Classes.OPdq import OPdq
from pyleecan.Classes.ForceMT import ForceMT
from pyleecan.Functions.Plot import dict_2D
from pyleecan.Classes.MagFEMM import MagFEMM
from os.path import join
from Tests import save_validation_path as save_path

# simple tests for geometry discretization
file_path = "C:/Users/pc/Downloads/easy_test_MagNetwork.json"
# file_path = "C:/Users/pc/Downloads/easy_test_MagNetwork2.json"
# file_path = "C:/Users/pc/Downloads/easy_test_MagNetwork3.json"

# Defining the validation cases paths
# file_path = "C:/Users/pc/Downloads/SPMSM_val.json"  # case study 1
# file_path = (
#     "C:/Users/pc/AppData/Roaming/pyleecan/Machine/Benchmark.json"  # case study 2
# )
# file_path = "C:/Users/pc/AppData/Roaming/pyleecan/Machine/SPMSM_18s16p_loss.json"  # case study 3
# file_path = (
#     "C:/Users/pc/AppData/Roaming/pyleecan/Machine/SPMSM_skew.json"  # case study 4
# )


# Loading the machine from the file path
SPMSM_val = load(file_path)

# Creating the simulation object of the MagNetwork
simu = Simu1(name="test_magnetwork", machine=SPMSM_val)

# simu.machine.stator.Rext - simu.machine.rotor.Rint

# Defining the Inputs and discretizations
# MagNetwork simulation
simu.input = InputCurrent(
    OP=OPdq(N0=1000, Id_ref=0, Iq_ref=0),
    Na_tot=90 * 4,  # simu.input.Na_tot
    Nt_tot=1,
)

# Defining the magnetic simulation: with periodicity for the MagNetwork
simu.mag = MagNetwork(
    is_periodicity_a=True,
    is_periodicity_t=False,
    type_model=1,
    type_coord_sys=2,
    Kmesh_fineness=2,
    rotor_shift=8,
)

# Cartesian meshing
# print(simu.mag.cartesianmeshclass_pyleecan())

# Test the MagNetwork simulation

# Initializing the harmonics calculation for simu
simu.force = ForceMT()

# Creating the simulation object of the FEMM
simu2 = Simu1(name="test_FEMM", machine=SPMSM_val)

# FEMM simulation
simu2.input = InputCurrent(
    OP=OPdq(N0=1000, Id_ref=0, Iq_ref=0),
    Na_tot=90 * 4,
    Nt_tot=1,
    is_periodicity_t=False,
    is_periodicity_a=True,
    angle_rotor_initial=1.0e-12,
)

# Defining the magnetic simulation: with periodicity for theFEMM
simu2.mag = MagFEMM(
    is_periodicity_a=True,
    is_periodicity_t=False,
    nb_worker=4,  # number of FEMM windows to be opened
    is_get_meshsolution=True,  # Get the mesh solution
    is_fast_draw=True,
    is_calc_torque_energy=False,
    Kmesh_fineness=2,  # Define the mesh fineness
)

# Initializing the harmonics calculation for simu2
simu2.force = ForceMT()

# Running both simulations
out = simu.run()  # MagNetwork
out2 = simu2.run()  # FEMM

# plotting the radial and tangential 2D flux density curve
# out.mag.B.plot_2D_Data("angle")
# out2.mag.B.plot_2D_Data("angle")

# Computing the harmonics
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

# # Plotting the results of the MagNetwork harmonics
# out.mag.B.plot_2D_Data(
#     # "time",
#     "angle[0]{°}",
#     # data_list=[out2.mag.B],
#     data_list=[out.mag.B],
#     legend_list=["Periodic", "Full"],
#     is_show_fig=True,
#     **dict_2D
# )

# computing the harmonics of the FEMM
# out2.mag.B.plot_2D_Data(
#     # "time",
#     "angle[0]{°}",
#     # data_list=[out2.mag.B],
#     data_list=[out2.mag.B],
#     legend_list=["Periodic", "Full"],
#     save_path=join(save_path, simu2.name + "_B_space_FEMM.png"),
#     is_show_fig=False,
#     **dict_2D
# )

# # Computing the harmonics using the Maxwell tensor of the MagNetwork class
# out.force.AGSF.plot_2D_Data(
#     "wavenumber=[0,10]",
#     "time[0]",
#     # data_list=[out2.force.AGSF],
#     data_list=[out.force.AGSF],
#     legend_list=["Periodic", "Full"],
#     is_show_fig=True,
#     **dict_2D
# )

# Computing the harmonics of the FEMM
# out2.force.AGSF.plot_2D_Data(
#     "wavenumber=[0,10]",
#     "time[0]",
#     data_list=[out2.force.AGSF],
#     legend_list=["Periodic", "Full"],
#     save_path=join(save_path, simu.name + "_P_space_fft_FEMM.png"),
#     is_show_fig=False,
#     **dict_2D
# )

# out.force.AGSF.plot_2D_Data(
#     "freqs",
#     "angle[0]",
#     data_list=[out2.force.AGSF],
#     legend_list=["Periodic", "Full"],
#     save_path=join(save_path, simu.name + "_B_space.png"),
#     is_show_fig=False,
#     **dict_2D
# )

# Comparison of the radial 2D flux density curve between MagNetwork and FEMM
out.mag.B.plot_2D_Data(
    "angle[smallestperiod]{°}",
    component_list=["radial"],
    data_list=[out2.mag.B],
    legend_list=["MagNetwork", "FEMM"],
)

# Comparison of the tangential 2D flux density curve between MagNetwork and FEMM
out.mag.B.plot_2D_Data(
    "angle[smallestperiod]{°}",
    component_list=["tangential"],
    data_list=[out2.mag.B],
    legend_list=["MagNetwork", "FEMM"],
)

# Comparison of the fft between MagNetwork and FEMM
out.mag.B.plot_2D_Data(
    "wavenumber=[0,10]",
    component_list=["radial"],
    data_list=[out2.mag.B],
    legend_list=["MagNetwork", "FEMM"],
)
pass
