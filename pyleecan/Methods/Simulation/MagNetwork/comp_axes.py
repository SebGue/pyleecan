from pyleecan.Classes.Magnetics import Magnetics


def comp_axes(self, output):

    # Call Magnetics.comp_axes()
    axes_dict = Magnetics.comp_axes(self, output)

    theta = np.linspace(
        axes_dict["angle"].initial,
        axes_dict["angle"].final,
        axes_dict["angle"].number + 1,
        # endpoint=False,
    )

    axes_dict["angle"] = (theta[1:] + theta[:-1]) / 2

    return axes_dict
