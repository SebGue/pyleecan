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
    mesh_pv=None,
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

    mesh_point = mesh_pv.cell_data_to_point_data()
    contours = mesh_point.contour()
    p.add_mesh(contours, color="black", line_width=3)