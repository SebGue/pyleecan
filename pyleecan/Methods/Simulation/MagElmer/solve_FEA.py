import numpy as np
import subprocess

from numpy import zeros, pi, floor_divide, loadtxt, sign, angle as np_angle
from os.path import join

from ....Functions.Winding.gen_phase_list import gen_name

from ....Functions.labels import (
    AIRGAP_LAB,
    ROTOR_LAB,
    ROTOR_LAB_S,
    short_label,
    decode_label,
    get_obj_from_label,
    LAM_LAB_S,
    STATOR_LAB_S,
    HOLEV_LAB_S,
    HOLEM_LAB_S,
    WIND_LAB_S,
    MAG_LAB,
    SHAFT_LAB,
    NO_LAM_LAB,
    SLID_LAB,
    BOT_LAB,
)
from ....Functions.Winding.find_wind_phase_color import get_phase_id
from .... import __version__

from ....Functions.get_path_binary import get_path_binary

from ....Classes.MachineSIPMSM import MachineSIPMSM
from ....Classes.MachineIPMSM import MachineIPMSM
from ....Methods import NotImplementedYetError


def solve_FEA(self, output, sym, angle, time, elmer_sif_file):
    """
    Solve Elmer model to calculate airgap flux density, torque instantaneous/average/ripple values,
    flux induced in stator windings and flux density, field and permeability maps

    Parameters
    ----------
    self: MagElmer
        A MagElmer object
    output: Output
        An Output object
    sym: int
        Spatial symmetry factor
    time: ndarray
        Time vector for calculation
    angle: ndarray
        Angle vector for calculation
    Is : ndarray
        Stator current matrix (qs,Nt) [A]
    Ir : ndarray
        Stator current matrix (qs,Nt) [A]
    angle_rotor: ndarray
        Rotor angular position vector (Nt,)
    elmer_sim_file: str
        Elmer solver input file
    """

    elmermesh_folder = self.get_path_save_fea(output)
    time = np.append(time, time[1] + time[-1])

    # setup Elmer solver
    # ElmerSolver v8.4 must be installed and in the PATH
    self.get_logger().debug("Solving Simulation")

    ElmerSolver_binary = get_path_binary("ElmerSolver")
    cmd_elmersolver = [ElmerSolver_binary, elmer_sif_file]
    self.get_logger().info(
        "Calling ElmerSolver: " + " ".join(map(str, cmd_elmersolver))
    )
    elmersolver = subprocess.Popen(
        cmd_elmersolver, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    (stdout, stderr) = elmersolver.communicate()
    elmersolver.wait()
    self.get_logger().info(stdout.decode("UTF-8"))
    if elmersolver.returncode != 0:
        self.get_logger().info("ElmerSolver [Error]: " + stderr.decode("UTF-8"))
        return False
    elmersolver.terminate()
    self.get_logger().info("ElmerSolver call complete!")

    self.get_meshsolution(output)

    Na = angle.size
    Nt = time.size - 1

    # Loading parameters for readibility
    L1 = output.simu.machine.stator.comp_length()

    scalars_file = join(elmermesh_folder, "scalars.dat")
    ecp, mfe, agt, iv, im, tq = loadtxt(
        scalars_file, unpack=True, usecols=(0, 1, 2, 3, 4, 5)
    )
    # ecp: eddy current power
    # mfe: magnetic field energy
    # agt: air gap torque
    # iv: inertial volume
    # im: inertial moment
    # tq: group 1 torque

    # TODO Load Air gap flux density

    # FEM_dict = output.mag.FEM_dict
    #
    if (
        hasattr(output.simu.machine.stator, "winding")
        and output.simu.machine.stator.winding is not None
    ):
        qs = output.simu.machine.stator.winding.qs  # Winding phase number
        Phi_wind_stator = zeros((Nt, qs))
    else:
        Phi_wind_stator = None

    # Initialize results matrix
    Br = zeros((Nt, Na))
    Bt = zeros((Nt, Na))
    Bz = zeros((Nt, Na))
    Tem = tq * sym  # torque per meter

    return Br, Bt, Bz, Tem, Phi_wind_stator
