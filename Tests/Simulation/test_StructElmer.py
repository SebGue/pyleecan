from os.path import join
from os import makedirs
from numpy import pi
import pytest
from Tests import save_validation_path as save_path, TEST_DATA_DIR
from pyleecan.Classes.OPdq import OPdq
from pyleecan.definitions import DATA_DIR

from pyleecan.Classes.LamSlotMag import LamSlotMag
from pyleecan.Classes.SlotM11 import SlotM11
from pyleecan.Classes.MachineSIPMSM import MachineSIPMSM
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.StructElmer import StructElmer
from pyleecan.Classes.InputVoltage import InputVoltage
from pyleecan.Classes.Output import Output
from pyleecan.Functions.load import load


# get the machine
machine_1 = load(join(DATA_DIR, "Machine", "Toyota_Prius.json"))
save_path = join(save_path, "StructElmer")

# mesh settings, original line label names have to be used (not the translated)
mm = 1e-3

mesh_dict_1 = {
    # define surface element sizes ...
    "Rotor-0_HoleMag_R0-T0-S0": 1 * mm,
    "Rotor-0_Lamination": 1 * mm,
    "Rotor-0_HoleVoid_R0-T2-S0": 0.5 * mm,
    # ... or number of elements per line (works only for lines with 'Boundary' property)
    "Rotor-0_YokeSide-Right-0": 100,
    "Rotor-0_LaminationYoke": 100,
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

        machine.rotor.slot = SlotM11(H0=0, W0=pi / 16, H1=1 * mm, W1=5 * mm)
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
    out = obj.test_StructElmer_HoleM50()
    out = obj.test_StructElmer_HoleM50_no_magnets()

    # test centrifugal force on a disc
    out = obj.test_StructElmer_disk()

    # # plot some results
    # out.struct.meshsolution.plot_deflection(label="disp", factor=20)
    # out.struct.meshsolution.plot_contour(label="disp")
    # out.struct.meshsolution.plot_mesh()
    print("Done")
