# -*- coding: utf-8 -*-
from os.path import isfile, join
import pytest

from pyleecan.Classes.HoleUD import HoleUD
from pyleecan.Classes.LamHole import LamHole
from pyleecan.Classes.SurfLine import SurfLine
from pyleecan.definitions import DATA_DIR
from pyleecan.Functions.load import load
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.MagFEMM import MagFEMM
from pyleecan.Classes.InputCurrent import InputCurrent
from pyleecan.Classes.Magnet import Magnet
from Tests import save_validation_path as save_path
import matplotlib.pyplot as plt

# For AlmostEqual
DELTA = 1e-4

mach_path = join(DATA_DIR, "Machine", "IPMSM_A.json")
if isfile(mach_path):
    IPMSM_A = load(mach_path)

    surf_list = IPMSM_A.rotor.hole[0].build_geometry()
    magnet_dict = {
        "magnet_0": IPMSM_A.rotor.hole[0].magnet_0,
        "magnet_1": IPMSM_A.rotor.hole[0].magnet_1,
    }

    HUD = HoleUD(Zh=4, surf_list=surf_list, magnet_dict=magnet_dict)

    IPMSM_B = IPMSM_A.copy()
    IPMSM_B.rotor.hole[0] = HUD

    surf_list = IPMSM_A.rotor.hole[0].build_geometry()
    magnet_dict = {
        "magnet_0": IPMSM_A.rotor.hole[0].magnet_0,
        "magnet_1": IPMSM_A.rotor.hole[0].magnet_1,
    }
    magnet_dict["magnet_0"].type_magnetization = 0
    magnet_dict["magnet_1"].type_magnetization = 0
    HUD2 = HoleUD(Zh=4, surf_list=surf_list, magnet_dict=magnet_dict)

    IPMSM_C = IPMSM_A.copy()
    IPMSM_C.rotor.hole[0] = HUD2
    IPMSM_C.rotor.is_stator = True


@pytest.mark.METHODS
class Test_HoleUD_meth(object):
    def test_comp_magnet_surface(self):
        """Check that the computation of the magnet surface"""
        exp = IPMSM_A.rotor.hole[0].comp_surface_magnets()
        result = IPMSM_B.rotor.hole[0].comp_surface_magnets()

        assert exp == pytest.approx(result, rel=0.01)

    def test_comp_surface(self):
        """Check that the computation of the slot surface is correct"""
        exp = IPMSM_A.rotor.hole[0].comp_surface()
        result = IPMSM_B.rotor.hole[0].comp_surface()

        assert exp == pytest.approx(result, rel=0.01)

    def test_build_geometry_mag(self):
        """check that curve_list is correct and Parallel"""
        assert IPMSM_C.rotor.hole[0].magnet_dict["magnet_0"].type_magnetization == 0
        assert IPMSM_C.rotor.hole[0].magnet_dict["magnet_1"].type_magnetization == 0
        IPMSM_C.rotor.hole[0].magnet_dict["magnet_0"].type_magnetization = 1
        IPMSM_C.rotor.hole[0].magnet_dict["magnet_1"].type_magnetization = 1

        surf_list = IPMSM_C.rotor.hole[0].build_geometry()
        assert len(surf_list) == 5
        assert surf_list[0].label == "Hole_Stator_R0_T0_S0"
        assert surf_list[1].label == "HoleMagnet_Stator_Parallel_N_R0_T0_S0"
        assert surf_list[2].label == "Hole_Stator_R0_T1_S0"
        assert surf_list[3].label == "HoleMagnet_Stator_Parallel_N_R0_T1_S0"
        assert surf_list[4].label == "Hole_Stator_R0_T2_S0"

        IPMSM_C.rotor.hole[0].magnet_dict["magnet_0"].type_magnetization = 0
        IPMSM_C.rotor.hole[0].magnet_dict["magnet_1"].type_magnetization = 0

    def test_build_geometry_no_mag(self):
        """check that curve_list is correct (Remove magnet)"""
        assert IPMSM_B.rotor.hole[0].magnet_dict["magnet_0"] is not None
        assert IPMSM_B.rotor.hole[0].magnet_dict["magnet_1"] is not None
        IPMSM_B.rotor.hole[0].remove_magnet()
        assert IPMSM_B.rotor.hole[0].magnet_dict["magnet_0"] is None
        assert IPMSM_B.rotor.hole[0].magnet_dict["magnet_1"] is None

        surf_list = IPMSM_B.rotor.hole[0].build_geometry()

        assert len(surf_list) == 5
        for ii, surf in enumerate(surf_list):
            assert type(surf) is SurfLine
            assert surf.label == "Hole_Rotor_R0_T" + str(ii) + "_S0"

        surf_list = IPMSM_C.rotor.hole[0].build_geometry()
        assert len(surf_list) == 5
        assert surf_list[0].label == "Hole_Stator_R0_T0_S0"
        assert surf_list[1].label == "HoleMagnet_Stator_Radial_N_R0_T0_S0"
        assert surf_list[2].label == "Hole_Stator_R0_T1_S0"
        assert surf_list[3].label == "HoleMagnet_Stator_Radial_N_R0_T1_S0"
        assert surf_list[4].label == "Hole_Stator_R0_T2_S0"

    def test_comp_surface_magnet_id(self):
        """check that ids are correct (Remove magnet)"""
        assert IPMSM_B.rotor.hole[0].comp_surface_magnet_id(0) == 0
        surf_list = IPMSM_C.rotor.hole[0].surf_list
        for surf in surf_list:
            surf.label = ""
        assert IPMSM_C.rotor.hole[0].comp_surface_magnet_id(0) == 0

    @pytest.mark.long
    @pytest.mark.MagFEMM
    def test_convert_UD(self):
        """Check that the simulation is the same with the usual hole and the UD equivalent"""
        IPMSM_A = load(mach_path)
        simu = Simu1(name="Simu_Hole_normal", machine=IPMSM_A)

        # Definition of the enforced output of the electrical module
        N0 = 2500
        Nt_tot = 1
        Na_tot = 2048

        simu.input = InputCurrent(
            Is=None,
            Ir=None,  # No winding on the rotor
            N0=N0,
            Nt_tot=Nt_tot,
            Na_tot=Na_tot,
            Id_ref=50,
            Iq_ref=0,
        )

        # Definition of the magnetic simulation (no symmetry)
        simu.mag = MagFEMM(
            type_BH_stator=0,
            type_BH_rotor=0,
            is_periodicity_a=True,
            Kgeo_fineness=0.75,
        )
        simu.force = None
        simu.struct = None

        simu_UD = simu.copy()
        hole = IPMSM_A.rotor.hole[0].convert_to_UD()
        assert len(hole.surf_list) == 5
        assert len(hole.magnet_dict) == 2
        assert "magnet_0" in hole.magnet_dict
        assert isinstance(hole.magnet_dict["magnet_0"], Magnet)
        assert "magnet_1" in hole.magnet_dict
        assert isinstance(hole.magnet_dict["magnet_1"], Magnet)
        assert hole.Zh == 8
        assert len(hole.magnetization_dict_offset) == 2
        assert "magnet_0" in hole.magnetization_dict_offset
        assert isinstance(hole.magnetization_dict_offset["magnet_0"], float)
        assert "magnet_1" in hole.magnetization_dict_offset
        assert isinstance(hole.magnetization_dict_offset["magnet_1"], float)
        simu_UD.machine.rotor.hole[0] = hole

        (R1, R2) = simu.machine.rotor.hole[0].comp_radius()
        (R3, R4) = simu_UD.machine.rotor.hole[0].comp_radius()
        assert R1 == pytest.approx(R3)
        assert R2 == pytest.approx(R4)
        out = simu.run()
        out2 = simu_UD.run()

        # Plot the result by comparing the two simulation
        plt.close("all")

        out.plot_2D_Data(
            "mag.B",
            "angle",
            data_list=[out2.mag.B],
            legend_list=["Normal", "User-Defined"],
            save_path=join(save_path, "test_HoleUD.png"),
            is_show_fig=False,
        )


if __name__ == "__main__":
    a = Test_HoleUD_meth()
    a.test_convert_UD()
