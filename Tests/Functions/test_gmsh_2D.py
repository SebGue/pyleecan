import pytest
import sys
import json

from os import makedirs, listdir
from os.path import join, isdir, split, splitext
from pyleecan.Functions.load import load
from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.Output import Output
from pyleecan.Classes.SlotM10 import SlotM10
from pyleecan.definitions import DATA_DIR
from Tests import save_plot_path
from pyleecan.Methods.Simulation.MagElmer import (
    MagElmer_BP_dict,
)

try:
    from pyleecan.Functions.GMSH.draw_GMSH import draw_GMSH
except:
    draw_GMSH = ImportError

mm = 1e-3

user_mesh_dict = {
    # surfaces: element size
    "None_Shaft": 4 * mm,
    "None_Airbox": 4 * mm,
    "Stator-0_SlidingBand": 0.5 * mm,
    "Rotor-0_SlidingBand": 0.5 * mm,
    # lines: number of elements
    "ABSide-Right": 6,
    "ABSide-Left": 6,
    "Stator-0_LaminationYoke": 150,
    "sliding_sideline_Top_Right": 2,
    "sliding_sideline_Top_Left": 2,
    "sliding_radius_Top": 100,
    "airgap_radius_Bot": 200,
    "sliding_sideline_Bot_Right": 4,
    "sliding_sideline_Bot_Left": 4,
    "sliding_radius_Bot": 200,
    "airgap_radius_Bot": 200,
}


def get_test_files(directory):
    return [
        join(directory, file) for file in listdir(directory) if file.endswith(".json")
    ]


@pytest.mark.long_5s
@pytest.mark.GMSH2D
@pytest.mark.parametrize("file_path", get_test_files(join(DATA_DIR, "Machine")))
def test_gmsh_machines(file_path):
    """Check generation of the 2D mesh with gmsh"""
    if isinstance(draw_GMSH, ImportError):
        raise ImportError("Fail to import draw_GMSH (gmsh package missing)")

    # Import the machine from a script
    machine = load(file_path)
    save_path = join(save_plot_path, "GMSH")
    if not isdir(save_path):
        makedirs(save_path)
    # Plot the machine
    # im = machine.plot()

    path, file = split(file_path)
    name = splitext(file)[0]

    # Create the Simulation
    mySimu = Simu1(name=name, machine=machine)
    myResults = Output(simu=mySimu)

    gmsh_dict = draw_GMSH(
        output=myResults,
        sym=machine.comp_periodicity_spatial()[0],
        boundary_prop=MagElmer_BP_dict,
        is_lam_only_S=False,
        is_lam_only_R=False,
        user_mesh_dict=user_mesh_dict,
        is_sliding_band=True,
        is_airbox=True,
        path_save=join(save_path, name + ".msh"),
    )

    with open(join(save_path, name + ".json"), "w") as fw:
        json.dump(gmsh_dict, fw, default=encode_complex, indent=4)

    return None


# @pytest.mark.long_5s
# # @pytest.mark.GMSH2D
# @pytest.mark.IPMSM
# def test_gmsh_ipm():
#     """Check generation of the 2D mesh with gmsh"""
#     if isinstance(draw_GMSH, ImportError):
#         raise ImportError("Fail to import draw_GMSH (gmsh package missing)")

#     # Import the machine from a script
#     Toyota_Prius = load(join(DATA_DIR, "Machine", "Toyota_Prius.json"))
#     Toyota_Prius.stator.slot.H1 = 1e-3
#     save_path = join(save_plot_path, "GMSH")
#     if not isdir(save_path):
#         makedirs(save_path)
#     # Plot the machine
#     # im = Toyota_Prius.plot()

#     # Create the Simulation
#     mySimu = Simu1(name="test_gmsh_ipm", machine=Toyota_Prius)
#     myResults = Output(simu=mySimu)

#     gmsh_dict = draw_GMSH(
#         output=myResults,
#         sym=8,
#         boundary_prop=MagElmer_BP_dict,
#         is_lam_only_S=False,
#         is_lam_only_R=False,
#         user_mesh_dict=user_mesh_dict,
#         is_sliding_band=True,
#         is_airbox=True,
#         path_save=join(save_path, "test_gmsh_ipm.msh"),
#     )

#     with open(join(save_path, "test_gmsh_ipm.json"), "w") as fw:
#         json.dump(gmsh_dict, fw, default=encode_complex, indent=4)

#     return None


# @pytest.mark.long_5s
# # @pytest.mark.GMSH2D
# @pytest.mark.SPMSM
# def test_gmsh_spm():
#     """Check generation of the 2D mesh with gmsh"""
#     if isinstance(draw_GMSH, ImportError):
#         raise ImportError("Fail to import draw_GMSH (gmsh package missing)")

#     # Import the machine from a script
#     PMSM_A = load(join(DATA_DIR, "Machine", "SPMSM_001.json"))
#     PMSM_A.rotor.slot = SlotM10(W1=15e-3, H1=3e-3, H0=0.0, W0=15e-3, Zs=8)

#     # PMSM_A.plot()
#     save_path = join(save_plot_path, "GMSH")
#     if not isdir(save_path):
#         makedirs(save_path)

#     # Create the Simulation
#     mySimu = Simu1(name="test_gmsh_spm", machine=PMSM_A)
#     myResults = Output(simu=mySimu)

#     gmsh_dict = draw_GMSH(
#         output=myResults,
#         sym=4,
#         boundary_prop=MagElmer_BP_dict,
#         is_lam_only_S=False,
#         is_lam_only_R=False,
#         user_mesh_dict=user_mesh_dict,
#         is_sliding_band=True,
#         is_airbox=True,
#         path_save=join(save_path, "test_gmsh_spm.msh"),
#     )

#     with open(join(save_path, "test_gmsh_spm.json"), "w") as fw:
#         json.dump(gmsh_dict, fw, default=encode_complex, indent=4)

#     return None


# @pytest.mark.long_5s
# # @pytest.mark.GMSH2D
# @pytest.mark.SPMSM
# def test_gmsh_benchmark():
#     """Check generation of the 2D mesh with gmsh"""
#     if isinstance(draw_GMSH, ImportError):
#         raise ImportError("Fail to import draw_GMSH (gmsh package missing)")

#     # Import the machine from a script
#     Benchmark = load(join(DATA_DIR, "Machine", "Benchmark.json"))
#     # Benchmark.stator.slot.H1 = 1e-3
#     save_path = join(save_plot_path, "GMSH")
#     if not isdir(save_path):
#         makedirs(save_path)
#     # Plot the machine
#     # im = Toyota_Prius.plot()

#     # Create the Simulation
#     mySimu = Simu1(name="test_gmsh_benchmark", machine=Benchmark)
#     myResults = Output(simu=mySimu)

#     gmsh_dict = draw_GMSH(
#         output=myResults,
#         sym=1,
#         boundary_prop=MagElmer_BP_dict,
#         is_lam_only_S=False,
#         is_lam_only_R=False,
#         user_mesh_dict=user_mesh_dict,
#         is_sliding_band=True,
#         is_airbox=True,
#         path_save=join(save_path, "test_gmsh_benchmark.geo"),
#     )

#     with open(join(save_path, "test_gmsh_benchmark.json"), "w") as fw:
#         json.dump(gmsh_dict, fw, default=encode_complex, indent=4)

#     return None


# @pytest.mark.long_5s
# # @pytest.mark.GMSH2D
# def test_gmsh_WRSM():
#     """Check generation of the 2D mesh with gmsh"""
#     if isinstance(draw_GMSH, ImportError):
#         raise ImportError("Fail to import draw_GMSH (gmsh package missing)")

#     # Import the machine from a script
#     Zoe = load(join(DATA_DIR, "Machine", "Renault_Zoe.json"))
#     # Benchmark.stator.slot.H1 = 1e-3
#     save_path = join(save_plot_path, "GMSH")
#     if not isdir(save_path):
#         makedirs(save_path)
#     # Plot the machine
#     # im = Toyota_Prius.plot()

#     # Create the Simulation
#     mySimu = Simu1(name="test_gmsh_Zoe", machine=Zoe)
#     myResults = Output(simu=mySimu)

#     gmsh_dict = draw_GMSH(
#         output=myResults,
#         sym=1,
#         boundary_prop=MagElmer_BP_dict,
#         is_lam_only_S=False,
#         is_lam_only_R=False,
#         user_mesh_dict=user_mesh_dict,
#         is_sliding_band=True,
#         is_airbox=True,
#         path_save=join(save_path, "test_gmsh_Zoe.msh"),
#     )

#     with open(join(save_path, "test_gmsh_Zoe.json"), "w") as fw:
#         json.dump(gmsh_dict, fw, default=encode_complex, indent=4)

#     return None


# @pytest.mark.long_5s
# # @pytest.mark.GMSH2D
# def test_gmsh_SCIM():
#     """Check generation of the 2D mesh with gmsh"""
#     if isinstance(draw_GMSH, ImportError):
#         raise ImportError("Fail to import draw_GMSH (gmsh package missing)")

#     # Import the machine from a script
#     RT = load(join(DATA_DIR, "Machine", "Railway_Traction.json"))
#     # Benchmark.stator.slot.H1 = 1e-3
#     save_path = join(save_plot_path, "GMSH")
#     if not isdir(save_path):
#         makedirs(save_path)
#     # Plot the machine
#     # im = Toyota_Prius.plot()

#     # Create the Simulation
#     mySimu = Simu1(name="test_gmsh_RT", machine=RT)
#     myResults = Output(simu=mySimu)

#     gmsh_dict = draw_GMSH(
#         output=myResults,
#         sym=1,
#         boundary_prop=MagElmer_BP_dict,
#         is_lam_only_S=False,
#         is_lam_only_R=False,
#         user_mesh_dict=user_mesh_dict,
#         is_sliding_band=True,
#         is_airbox=True,
#         path_save=join(save_path, "test_gmsh_RT.msh"),
#     )


def encode_complex(z):
    if isinstance(z, complex):
        return (z.real, z.imag)


if __name__ == "__main__":
    # gmsh_dict = test_gmsh_ipm()
    # gmsh_dict = test_gmsh_spm()
    # test_gmsh_benchmark()
    # test_gmsh_SCIM()
    pass
