# -*- coding: utf-8 -*-
from os.path import join

from ....Functions.GMSH.draw_GMSH import draw_GMSH
from ....Classes.OutMagElmer import OutMagElmer
from ....Methods.Simulation.MagElmer import MagElmer_BP_dict


def comp_flux_airgap(self, output, axes_dict, Is_val=None, Ir_val=None):
    """Build and solve Elmer model to calculate and store magnetic quantities

    Parameters
    ----------
    self : MagElmer
        a MagElmer object
    output : Output
        an Output object
    axes_dict: {Data}
        Dict of axes used for magnetic calculation
    """

    # Init output dict
    out_dict = dict()
    if output.mag.internal is None:
        output.mag.internal = OutMagElmer()

    # Get time and angular axes
    Angle = axes_dict["angle"]
    Time = axes_dict["time"]

    # Set the angular symmetry factor according to the machine and check if it is anti-periodic
    sym, is_antiper_a = Angle.get_periodicity()

    # Import angular vector from Data object
    angle = Angle.get_values(
        is_oneperiod=self.is_periodicity_a,
        is_antiperiod=is_antiper_a and self.is_periodicity_a,
    )
    # Na = angle.size

    # Check if the time axis is anti-periodic
    _, is_antiper_t = Time.get_periodicity()

    # Number of time steps
    time = Time.get_values(
        is_oneperiod=self.is_periodicity_t,
        is_antiperiod=is_antiper_t and self.is_periodicity_t,
    )
    Nt = time.size

    # Get rotor angular position
    angle_rotor = output.get_angle_rotor()[0:Nt]

    # Setup the Elmer simulation
    # Geometry building
    gmsh_filename = self.get_path_save_fea(output) + ".msh"
    if not self.import_file:  # True if None or len == 0
        self.get_logger().debug("Drawing machine in GMSH...")
        output.mag.internal.FEA_dict = draw_GMSH(
            output=output,
            sym=sym,
            boundary_prop=MagElmer_BP_dict,
            is_lam_only_S=False,
            is_lam_only_R=False,
            user_mesh_dict=self.FEA_dict,
            is_sliding_band=True,
            is_airbox=True,
            path_save=gmsh_filename,
        )

    else:
        self.get_logger().debug("Reusing the FEA file: " + self.import_file)
        # output.mag.internal.FEA_dict = self.FEA_dict
        pass

    # post process GMSH mesh with ElmerGrid
    if not self.gen_elmer_mesh(output):
        self.get_logger().error(
            "Something went wrong while processing mesh with ElmerGrid."
        )

    # generate the elmer solver input file
    elmer_sif_file = self.gen_elmer_sif(output, sym, time, angle_rotor, Is_val, Ir_val)

    # Solve for all time step and store all the results in output
    Br, Bt, Bz, Tem, Phi_wind_stator = self.solve_FEA(
        output, sym, angle, time, elmer_sif_file
    )

    # Store standards Magnetics outputs in out_dict
    out_dict["B_{rad}"] = Br
    out_dict["B_{circ}"] = Bt
    out_dict["B_{ax}"] = Bz
    out_dict["Tem"] = Tem
    out_dict["Phi_wind_stator"] = Phi_wind_stator

    # TODO store other Elmer outputs in out_dict

    return out_dict
