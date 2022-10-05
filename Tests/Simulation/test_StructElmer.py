from PySide2 import QtWidgets, QtGui, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvas  # for proper plots

from os.path import join
from os import makedirs
from numpy import pi
import pytest
from Tests import save_validation_path as save_path, TEST_DATA_DIR
from pyleecan.Classes.OPdq import OPdq
from pyleecan.definitions import DATA_DIR

from pyleecan.Classes.LamSlotMag import LamSlotMag
from pyleecan.Classes.SlotM11 import SlotM11
from pyleecan.Classes.HoleM52R import HoleM52R
from pyleecan.Classes.BoreSinePole import BoreSinePole
from pyleecan.Classes.MachineSIPMSM import MachineSIPMSM
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.StructElmer import StructElmer
from pyleecan.Classes.InputVoltage import InputVoltage
from pyleecan.Classes.Output import Output
from pyleecan.Functions.load import load
from pyleecan.Functions.MeshSolution.get_indices_limited import get_indices_limited
from pyleecan.Functions.MeshSolution.get_area import get_area

# get the machine
machine_1 = load(join(DATA_DIR, "Machine", "Toyota_Prius.json"))
save_path = join(save_path, "StructElmer")
# mesh settings, original line label names have to be used (not the translated)
n1 = 6
n2 = 40

mesh_dict_1 = {
    "Rotor": {
        "Rotor_Magnet_Top_0": n2,
        "Rotor_Magnet_Bottom_0": n2,
        "Rotor_Magnet_Left_0": n1,
        "Rotor_Magnet_Right_0": n1,
        "Rotor_Magnet_Top_1": n2,
        "Rotor_Magnet_Bottom_1": n2,
        "Rotor_Magnet_Left_1": n1,
        "Rotor_Magnet_Right_1": n1,
        "Rotor_Hole_Top_0": n2,
        "Rotor_Hole_Left_0": n1,
        "Rotor_Hole_Right_0": n1,
        "Rotor_Hole_Top_1": 0,
        "Rotor_Hole_Left_1": n1,
        "Rotor_Hole_Right_1": n1,
        "Rotor_Tangential_Bridge": 40,
        "Rotor_Radial_Bridge": 40,
        "ROTOR_BORE_CURVE": 100,
        "Lamination_Rotor_Bore_Radius_Ext": 100,
    }
}


@pytest.mark.long_5s
@pytest.mark.StructElmer
@pytest.mark.IPMSM
@pytest.mark.SingleOP
class Test_StructElmer(object):
    """Test some basic workflow of StructElmer simulations"""

    def test_StructElmer_HoleM50(self):
        """Test StructElmer simulation with 2 magnets on HoleM50 rotor"""

        # copy the machine
        machine = machine_1.copy()

        # some modifications to geometry
        machine.rotor.hole[0].W2 = 1.0e-3

        # setup the simulation
        simu = Simu1(name="test_StructElmer_HoleM50", machine=machine)
        output = Output(simu=simu)
        output.path_result = join(save_path, "Hole50")
        makedirs(output.path_result)

        simu.struct = StructElmer()
        simu.struct.FEA_dict_enforced = mesh_dict_1
        simu.struct.is_get_mesh = True

        # set rotor speed and run simulation
        simu.input = InputVoltage(OP=OPdq(N0=10000))  # rpm
        simu.run()

        return output

    def test_StructElmer_HoleM52R(self):
        """Test StructElmer simulation with 1 magnet on HoleM52R rotor"""
        mm = 1e-3

        # copy the machine
        machine = machine_1.copy()
        machine.name = "Prius stator - HoleM52R rotor"

        # some modifications to geometry
        hole = HoleM52R()
        bore = BoreSinePole()
        machine.rotor.hole[0] = hole
        machine.rotor.bore = bore

        hole.Zh = 8
        hole.W0 = 40 * mm
        hole.W1 = 2 * mm
        hole.R0 = 1 * mm
        hole.H0 = 2 * mm
        hole.H1 = 5 * mm
        hole.H2 = 1 * mm

        bore.W0 = 46 * mm
        bore.delta_d = 0.75 * mm
        bore.delta_q = 8.00 * mm

        machine.plot()
        # setup the simulation
        simu = Simu1(name="test_StructElmer_HoleM52R", machine=machine)
        output = Output(simu=simu)
        output.path_result = join(save_path, "Hole52R")
        makedirs(output.path_result)

        simu.struct = StructElmer()
        simu.struct.FEA_dict_enforced = mesh_dict_1
        simu.struct.is_get_mesh = True

        # set rotor speed and run simulation
        simu.input = InputVoltage(OP=OPdq(N0=20000))  # rpm
        simu.run()

        return output

    def test_StructElmer_HoleM50_no_magnets(self):
        """Test StructElmer simulation without magnets on HoleM50 rotor"""

        # get the machine
        machine = machine_1.copy()

        # some modifications to geometry
        machine.rotor.hole[0].W2 = 1.0e-3
        # machine.rotor.hole[0].H2 = 0.0e-3
        # machine.rotor.hole[0].W1 = 1.0e-3

        # setup the simulation
        simu = Simu1(name="test_StructElmer_HoleM50_no_magnets", machine=machine)
        output = Output(simu=simu)
        output.path_result = join(save_path, "Hole50_no_mag")
        makedirs(output.path_result)

        simu.struct = StructElmer()
        simu.struct.FEA_dict_enforced = mesh_dict_1
        simu.struct.include_magnets = False
        simu.struct.is_get_mesh = True

        # set rotor speed and run simulation
        simu.input = InputVoltage(OP=OPdq(N0=10000))  # rpm
        simu.run()

        return output

    def test_StructElmer_disk(self):
        """Test StructElmer simulation with disc geometry (i.e. slotless rotor)"""
        # TODO compare to analytical values

        # setup new machine and copy stator props of ref. machine
        machine = MachineSIPMSM()
        machine.stator = machine_1.stator.copy()
        machine.rotor = LamSlotMag()

        machine.rotor.Rint = machine_1.rotor.Rint
        machine.rotor.Rext = machine_1.rotor.Rext
        machine.rotor.mat_type = machine_1.rotor.mat_type.copy()

        machine.rotor.slot = SlotM11(H0=0, W0=pi / 16)
        machine.rotor.slot.Zs = 8
        machine.rotor.is_stator = False

        # setup the simulation
        simu = Simu1(name="test_StructElmer_disk", machine=machine)
        output = Output(simu=simu)
        output.path_result = join(save_path, "disk")
        makedirs(output.path_result)

        simu.struct = StructElmer()
        # simu.struct.FEA_dict_enforced = mesh_dict_1
        simu.struct.include_magnets = False
        simu.struct.is_get_mesh = True

        # set rotor speed and run simulation
        simu.input = InputVoltage(OP=OPdq(N0=10000))  # rpm
        simu.run()

        return output


# To run it without pytest
if __name__ == "__main__":
    # create test object
    obj = Test_StructElmer()
    # test Toyota_Prius (HoleM50-Rotor) with minor modification
    # out = obj.test_StructElmer_HoleM50()
    # out = obj.test_StructElmer_HoleM50_no_magnets()

    out = obj.test_StructElmer_HoleM52R()

    # test centrifugal force on a disc
    # out = obj.test_StructElmer_disk()

    # get_indices_limited(out.struct.meshsolution,label="vonmises")
    field = out.struct.meshsolution.get_field(
        label="vonmises",
        index=None,
        indices=None,
        is_rthetaz=False,
        is_pol2cart=False,
        is_radial=False,
        is_normal=False,
        is_rms=False,
        is_center=False,
        is_surf=False,
        is_squeeze=True,
    )

    # # plot some results
    # out.struct.meshsolution.plot_deflection(label="disp", factor=20)
    out.struct.meshsolution.plot_contour(label="disp")
    out.struct.meshsolution.plot_contour(label="vonmises")
    out.struct.meshsolution.plot_mesh()
