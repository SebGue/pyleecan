import subprocess
import re

from matplotlib import pyplot

from numpy import zeros, loadtxt
from os.path import split
from os.path import join

from .... import __version__
from ....Functions.get_path_binary import get_path_binary
from ....Functions.labels import (
    STATOR_LAB_S,
    WIND_LAB_S,
    decode_label,
    get_obj_from_label,
)
from ....Functions.Winding.find_wind_phase_color import get_phase_id


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
    # readability
    machine = output.simu.machine

    # setup Elmer solver
    # ElmerSolver v8.4 must be installed and in the PATH
    # TODO enable parallel computation

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
    L1 = machine.stator.comp_length()
    save_path = self.get_path_save(output)

    scalars_file = join(elmermesh_folder, "scalars.dat")
    results = _get_scalars(scalars_file)

    # get data
    tq = results["group 1 torque"]
    agt = results["air gap torque"]

    # ecp: eddy current power
    # mfe: magnetic field energy
    # iv: inertial volume
    # im: inertial moment

    # get flux linkage
    if hasattr(machine.stator, "winding") and machine.stator.winding is not None:
        qs = machine.stator.winding.qs  # Winding phase number

        # get the stator winding data from results
        wind_dict = {}
        for label, item in results.items():
            if label.startswith(STATOR_LAB_S.lower()) and WIND_LAB_S.lower() in label:
                # fix lower case Elmer names
                label = label.replace(WIND_LAB_S.lower(), WIND_LAB_S)
                label = label.replace(STATOR_LAB_S.lower(), STATOR_LAB_S)
                # get and store data
                phase_id, Ncond = _get_wind_props(machine, label)
                wind_dict[label] = {}
                wind_dict[label]["results"] = item
                wind_dict[label]["phase"] = phase_id
                wind_dict[label]["Ncond"] = Ncond

        Phi_wind_stator = zeros((Nt, qs))
        for item in wind_dict.values():
            Phi_wind_stator[:, item["phase"]] += item["results"] * item["Ncond"]

    else:
        Phi_wind_stator = None

    pyplot.plot(time, Phi_wind_stator)
    pyplot.show()

    # get the air gap flux result
    Br = zeros((Nt, Na))
    Bt = zeros((Nt, Na))
    Bz = zeros((Nt, Na))

    # get the torque
    # TODO compare different Elmer torque calc. methods
    Tem = tq * sym * L1

    return Br, Bt, Bz, Tem, Phi_wind_stator


def _get_scalars(scalars_file):
    """
    Read the Elmer scalar result files and return all data as a dict with
    respective variable names as keys.
    """
    scalars_names_file = scalars_file + ".names"
    scalars = loadtxt(scalars_file, unpack=True)

    # load names
    with open(scalars_names_file, "r") as file:
        file_data = file.read()

    # extract variable names
    lines = file_data.splitlines()
    variable_lines = [line for line in lines if re.match(r"^\s*\d+:\s+", line)]

    # init variable names
    data = {}

    # process lines
    for line in variable_lines:
        col_number, typ, var = [txt.strip() for txt in line.split(":")]
        var = var.replace("a mask ", "", 1)
        data[var] = scalars[int(col_number) - 1, :]

    return data


def _get_wind_props(machine, label):
    """Get some winding properties by surface label."""
    label_dict = decode_label(label)
    lam_obj = get_obj_from_label(machine, label_dict=label_dict)
    wind_mat = lam_obj.winding.get_connection_mat(lam_obj.get_Zs())
    Nrad_id = label_dict["R_id"]  # zone radial coordinate
    Ntan_id = label_dict["T_id"]  # zone tangential coordinate
    Zs_id = label_dict["S_id"]  # Zone slot number coordinate
    # Get the phase value in the correct slot zone
    phase_id = get_phase_id(wind_mat, Nrad_id, Ntan_id, Zs_id)
    Ncond = wind_mat[Nrad_id, Ntan_id, Zs_id, phase_id]

    return phase_id, Ncond
