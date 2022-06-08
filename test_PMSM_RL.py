from os.path import join

import numpy as np
from pyleecan.Classes.InputCurrent import InputCurrent
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.OPdq import OPdq
from pyleecan.Classes.MagFEMM import MagFEMM


from pyleecan.Functions.load import load

from pyleecan.definitions import DATA_DIR

is_show_fig = True


def test_FEMM_SPMSM():

    machine = load(join(DATA_DIR, "Machine", "SPMSM_264s44p.json"))

    simu = Simu1(name="test_FEMM", machine=machine)

    # Ic = 230 * np.exp(1j * 140 * np.pi / 180)
    Ic = 0
    SPEED = 0

    simu.input = InputCurrent(
        Nt_tot=10 * 44,  # time steps
        Na_tot=2000 * 44,  # number of points in the air gap (angle discretization)
        OP=OPdq(N0=SPEED, Id_ref=Ic.real, Iq_ref=Ic.imag),  # operating point
        is_periodicity_t=True,
        is_periodicity_a=True,
    )

    simu.mag = MagFEMM(
        is_periodicity_a=True,
        is_periodicity_t=True,
        nb_worker=4,  # nbre de fenetre de feem ouvertes
        is_get_meshsolution=True,
        is_fast_draw=True,
        is_calc_torque_energy=False,
    )

    out = simu.run()

    if is_show_fig:
        # plot the air gap flux density on one line only
        out.mag.B.plot_2D_Data("angle", "time[1]", component_list=["radial"])
        # out.mag.B.plot_2D_Data("angle", "time[1]", component_list=["tangential"])

        #
        out.mag.meshsolution.plot_contour(
            label="B",
            is_show_fig=True,
        )

    return out


# To run it without pytest
if __name__ == "__main__":

    test_FEMM_SPMSM()
