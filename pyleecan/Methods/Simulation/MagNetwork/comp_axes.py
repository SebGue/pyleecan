from pyleecan.Classes.Magnetics import Magnetics
import numpy as np
from SciDataTool import Data1D, DataLinspace


def comp_axes(self, output):
    """Compute the axes required in the MagNetwork class

    Parameters
    ----------
    self : MagNetwork
        a MagNetwork object
    output : Output
        an Output object (to update)

    Returns
    -------
    axes_dict: {Data}
        Dict containing the dual-angle axis used in the MagNetwork class
    """

    # Call Magnetics.comp_axes()
    axes_dict = Magnetics.comp_axes(self, output)

    # Define the theta axis, taking into account the final point
    angle = axes_dict["angle"].get_values(is_smallestperiod=True)
    dtheta = angle[1] - angle[0]
    theta = np.linspace(
        axes_dict["angle"].initial,
        (axes_dict["angle"].final) + dtheta,
        axes_dict["angle"].number + 1,
    )

    # Save the original angle (theta_primal) of size N_point_theta
    axes_dict["theta_primal"] = Data1D(
        name="angle",
        unit="rad",
        values=theta,
        symmetries=axes_dict["angle"].symmetries.copy(),
    )

    # Definition of the dual angle
    angle_dual = (theta[1:] + theta[:-1]) / 2

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
