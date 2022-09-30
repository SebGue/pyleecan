from matplotlib import cm
from numpy import (
    real,
    min as np_min,
    max as np_max,
    abs as np_abs,
)
from ....definitions import config_dict
import pyvista as pv

# COLOR_MAP = config_dict["PLOT"]["COLOR_DICT"]["COLOR_MAP"]
COLOR_MAP = cm.get_cmap("hsv")


def plot_mesh_field(
    p,
    sargs,
    field_name,
    clim=None,
    mesh_pv=1,
    field=None,
    phase=1,
):

    mesh_pv[field_name] = real(field * phase)
    mesh_field = mesh_pv

    p.add_mesh(
        mesh_field,
        scalars=field_name,
        show_edges=False,
        cmap=COLOR_MAP,
        clim=clim,
        scalar_bar_args=sargs,
    )

    # streamlines = mesh_pv.streamlines(n_points=40, source_center=(2.35, 0.08, 0))
    # mesh_point = mesh_pv.cell_data_to_point_data()
    # streamlines = mesh_point.streamlines_evenly_spaced_2D(
    #     start_position=(2.35, 0.08, 0.0),
    #     separating_distance=3,
    #     separating_distance_ratio=0.2,
    #     compute_vorticity=False,  # vorticity already exists in dataset
    # )
    # p.add_mesh(streamlines.tube(radius=0.01), color="black")
