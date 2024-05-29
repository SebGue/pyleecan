import subprocess

from numpy import zeros, loadtxt
from os.path import split
from os.path import join

from .... import __version__
from ....Functions.get_path_binary import get_path_binary


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

    # setup Elmer solver
    # ElmerSolver v8.4 must be installed and in the PATH

    elmermesh_folder = self.get_path_save_fea(output)
    project_name = split(elmer_sif_file)[0]
    ElmerSolver_binary = get_path_binary("ElmerSolver")

    cmd_elmersolver = [ElmerSolver_binary, split(elmer_sif_file)[1]]

    self.get_logger().info(
        "Calling ElmerSolver: " + " ".join(map(str, cmd_elmersolver))
    )

    elmersolver = subprocess.Popen(
        cmd_elmersolver,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=project_name,
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
    Nt = time.size

    # Loading parameters for readibility
    L1 = output.simu.machine.stator.comp_length()
    save_path = self.get_path_save(output)

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
    Tem = tq * sym * L1

    # Phi_wind_stator = zeros((Nt, qs))

    # compute the data for each time step
    # TODO Other than FEMM, in Elmer I think it's possible to compute
    #      all time steps at once
    self.get_logger().debug("Solving Simulation")

    # run the computation
    if self.nb_worker > 1:
        # TODO run solver in parallel
        pass
    else:
        # TODO run solver 'normal'
        pass

    # get the air gap flux result
    # TODO add function (or method)
    # ii -> Time, jj -> Angle
    # Br[ii, jj], Bt[ii, jj] = get_airgap_flux()

    # get the torque
    # TODO add function (or method)
    # Tem[ii] = comp_Elmer_torque(FEM_dict, sym=sym)

    # flux linkage computation
    # if Phi_wind_stator is not None:
    #     # TODO
    #     # Phi_wind[ii, :] = comp_Elmer_Phi_wind()
    #     pass

    return Br, Bt, Bz, Tem, Phi_wind_stator
