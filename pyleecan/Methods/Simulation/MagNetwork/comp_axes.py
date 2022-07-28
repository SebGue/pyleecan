from pyleecan.Classes.Magnetics import Magnetics
import numpy as np
from SciDataTool import Data1D, DataLinspace


def comp_axes(self, output):

    # Call Magnetics.comp_axes()
    axes_dict = Magnetics.comp_axes(self, output)

    angle = axes_dict["angle"].get_values(is_smallestperiod=True)
    dtheta = angle[1] - angle[0]
    theta = np.linspace(
        axes_dict["angle"].initial,
        (axes_dict["angle"].final) + dtheta,
        axes_dict["angle"].number + 1,
        # endpoint=True,
    )

    # Save original angle (theta_primal) of size N_point_theta
    axes_dict["theta_primal"] = Data1D(
        name="angle",
        unit="rad",
        # values=angle_dual,
        values=theta,
        symmetries=axes_dict["angle"].symmetries.copy(),
    )

    angle_dual = (theta[1:] + theta[:-1]) / 2

    # axes_dict = {"angle": (theta[1:] + theta[:-1]) / 2}
    # axes_dict["angle"] = Data1D(
    #     name="angle",
    #     unit="rad",
    #     values=angle_dual,
    #     # values=theta,
    #     symmetries=axes_dict["angle"].symmetries.copy(),
    # )
    axes_dict["angle"] = DataLinspace(
        name="angle",
        unit="rad",
        initial=angle_dual[0],
        final=angle_dual[-1],
        number=len(angle_dual),
        include_endpoint=True,
        # values=theta,
        symmetries=axes_dict["angle"].symmetries.copy(),
    )
    return axes_dict
