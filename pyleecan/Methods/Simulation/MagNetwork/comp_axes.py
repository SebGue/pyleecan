from pyleecan.Classes.Magnetics import Magnetics


def comp_axes(self, output):

    # Call Magnetics.comp_axes()
    Magnetics.comp_axes(self, output)

    # theta ??

    angle = (theta[1:] + theta[:-1]) / 2

    axes_dict = output.simu.input.comp_axes(
        axes_list=["time", "angle"],
        axes_dict_in=axes_dict_geo,
        is_periodicity_a=self.is_periodicity_a,
        is_periodicity_t=self.is_periodicity_t,
        is_periodicity_rotor=self.is_periodicity_rotor,
    )

    return axes_dict
