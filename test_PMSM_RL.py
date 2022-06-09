from os.path import join

import numpy as np
from pyleecan.Classes.InputCurrent import InputCurrent
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.OPdq import OPdq
from pyleecan.Classes.MagFEMM import MagFEMM


from pyleecan.Functions.load import load

from pyleecan.definitions import DATA_DIR

from matplotlib import pyplot as plt
import numpy as np

# from qtpy.QtWidgets import *

is_show_fig = True


def test_FEMM_SPMSM():

    machine = load(join(DATA_DIR, "Machine", "SPMSM_264s44p.json"))

    # Import the simulation
    simu = Simu1(name="test_FEMM", machine=machine)

    # Ic = 230 * np.exp(1j * 140 * np.pi / 180)

    # Definition of the current ans speed
    Ic = 0
    SPEED = 0

    # Definition of time, angular discretization and operating point
    simu.input = InputCurrent(
        Nt_tot=10 * 44,  # time steps (time discretization)
        Na_tot=2000 * 44,  # number of points in the air gap (angle discretization)
        OP=OPdq(N0=SPEED, Id_ref=Ic.real, Iq_ref=Ic.imag),  # operating point
        # Periodicity
        is_periodicity_t=True,
        is_periodicity_a=True,
        angle_rotor_initial=0,
    )

    simu.mag = MagFEMM(
        is_periodicity_a=True,
        is_periodicity_t=True,
        nb_worker=4,  # number of FEMM windows to be opened
        is_get_meshsolution=True,  # Get the mesh solution
        is_fast_draw=True,
        is_calc_torque_energy=False,
        Kmesh_fineness=10,  # Define the mesh fineness
    )

    out = simu.run()

    # Get the radial B_rms
    Br_rms = out.mag.B.get_rphiz_along("time=rms", "angle=rms")["radial"]

    # Get the harmonics of a defined number of harmonics
    harmonics_number = 10
    result_harmonics = out.mag.B.components["radial"].get_harmonics(
        harmonics_number, "freqs", "wavenumber"
    )  # to extract 10 largest harmonics

    Br_harm = result_harmonics["B_{rad}"]
    freqs = result_harmonics["freqs"]
    wavenumber = result_harmonics["wavenumber"]

    # Get B for one defined harmonic
    f_specific = 50
    harmonic_specific = 3
    Br_specific = out.mag.B.components["radial"].get_magnitude_along(
        "freqs= 50/3", "wavenumber=3"
    )[
        "B_{rad}"
    ]  # to extract one specific harmonic

    harmonics = {
        "Bradial_harmonics": Br_harm,
        "harmonics_frequencies": freqs,
        "harmonics_wavenumber": wavenumber,
        "Br_specific_harmonic": Br_specific,
    }

    # Get electromagnetic torque
    Tem_whole = out.mag.Tem.get_along("time")[
        "T_{em}"
    ]  # reconstructed on whole time axis

    Tem_period = out.mag.Tem.get_along("time[smallestperiod]")[
        "T_{em}"
    ]  # on one period

    EM_torque = {"Tem_whole_periods": Tem_whole, "Tem_one_period": Tem_period}

    # Plot the B magnetic field mapping and the air gap flux density curve
    if is_show_fig:

        # plot the air gap flux density on one line only
        # out.mag.B.plot_2D_Data("angle", "time[1]", component_list=["radial"])

        # Plot the magnetic flux density mapping
        out.mag.meshsolution.plot_contour(
            # "time[1]",
            label="B",
            is_show_fig=True,
            # win_title=True,
            title="Flux Density Mapping B",
            win_title="Flux_Density_B",
        )

        # Plot the mesh of the machine under study
        out.mag.meshsolution.plot_mesh(
            # group_names=["stator core", "airgap", "rotor core"],
            is_show_fig=True,
        )

        # Plot the flux field mapping (flux density + flux lines)
        out.plot_B_mesh(
            is_show_fig=True,
            title="Flux Density and lines Mapping",
            win_title="Flux_Density_B",
        )

    return out, Br_rms, harmonics, EM_torque


# To run it without pytest
if __name__ == "__main__":

    test_FEMM_SPMSM()