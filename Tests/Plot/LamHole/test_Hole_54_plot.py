# -*- coding: utf-8 -*-

import pytest

from os.path import join
import matplotlib.pyplot as plt
from numpy import pi

from pyleecan.Classes.Frame import Frame
from pyleecan.Classes.LamHole import LamHole
from pyleecan.Classes.LamSlotWind import LamSlotWind
from pyleecan.Classes.MachineIPMSM import MachineIPMSM
from pyleecan.Classes.Magnet import Magnet
from pyleecan.Classes.Shaft import Shaft
from pyleecan.Classes.HoleM54 import HoleM54
from Tests import save_plot_path as save_path
from pyleecan.Functions.labels import BOUNDARY_PROP_LAB


class Test_Hole_54_plot(object):
    """pytest for Lamination with Hole plot"""

    @pytest.fixture
    def machine(self):
        """Run at the begining of every test to setup the machine"""
        plt.close("all")
        test_obj = MachineIPMSM()
        test_obj.rotor = LamHole(
            is_internal=True, Rint=0.1, Rext=0.2, is_stator=False, L1=0.7
        )
        test_obj.rotor.hole = list()
        test_obj.rotor.hole.append(
            HoleM54(Zh=8, W0=pi / 4, H0=50e-3, H1=10e-3, R1=100e-3)
        )
        test_obj.rotor.hole.append(
            HoleM54(Zh=8, W0=pi / 6, H0=25e-3, H1=10e-3, R1=100e-3)
        )

        return test_obj

    def test_Lam_Hole_54_plot(self, machine):
        """Test machine plot hole 54"""

        machine.rotor.plot(is_show_fig=False)
        fig = plt.gcf()
        fig.savefig(join(save_path, "test_Lam_Hole_s54-Rotor.png"))
        assert len(fig.axes[0].patches) == 18

        machine.rotor.hole[0].plot(is_show_fig=False)
        fig = plt.gcf()
        fig.savefig(join(save_path, "test_Lam_Hole_s54-Rotor hole.png"))
        assert len(fig.axes[0].patches) == 1

    def test_plot_line_labels(self, machine):
        """Test if the line labels are assigned to the respective lines."""
        machine.rotor.plot(is_show_fig=False)
        fig = plt.gcf()

        surfs = machine.rotor.build_geometry()

        for surf in surfs:
            lines = surf.get_lines()
            for line in lines:
                mid = line.get_middle()
                label = (
                    line.prop_dict.get(BOUNDARY_PROP_LAB, None)
                    if line.prop_dict
                    else None
                )
                if label:
                    plt.text(mid.real, mid.imag, label, fontsize=1)
        fig.savefig(join(save_path, "test_Lam_Hole_s54_line_label.png"), dpi=1000)

    def test_plot_point_labels(self, machine):
        """Test if the point labels are assigned to the respective points."""
        machine.rotor.hole[0].plot_schematics(
            is_add_point_label=True,
            is_add_schematics=False,
            is_add_main_line=False,
            is_show_fig=False,
        )
        fig = plt.gcf()
        fig.savefig(join(save_path, "test_Lam_Hole_s54_point_label.png"), dpi=1000)
