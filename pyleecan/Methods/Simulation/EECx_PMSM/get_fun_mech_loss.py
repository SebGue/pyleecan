# -*- coding: utf-8 -*-
from numpy import linspace, zeros, pi
import scipy.interpolate as interp


def get_fun_mech_loss(self, max_speed=50000):
    """Method to get the mech. losses (torque) function. The given mech. loss models
    are assumed to be independend of current. The sum of mech. losses will be used
    to create a loss interpolation function.
    """
    # some settings
    NN = 51  # number of speed interpolates

    # ref. speed vector for loss calculation
    speed_ref = linspace(0, max_speed, NN).tolist()

    if self.mech_loss is None or not self.mech_loss:
        return lambda Id, Iq, freq: 0

    losses = zeros(NN)
    for loss_mdl in self.mech_loss.values():
        # store original value for speed and set desired speed
        loss_mdl = loss_mdl.copy()
        loss_mdl.N0 = speed_ref
        loss_data = loss_mdl.comp_loss(output=self.table)[0]

        if loss_data is not None:
            # get data axes to compute means and sums over the axes
            axes = loss_data.get_axes()
            axes_list = [
                ax.name + "=sum"
                for ax in axes
                if ax.name not in ["time", "freqs", "speed"]
            ]

            data_dict = loss_data.get_along("time=mean", *axes_list, "speed")
            data = data_dict[loss_data.symbol]
            speed = data_dict["speed"]

        losses += data

    t_fric = zeros(NN)
    ii = speed != 0
    t_fric[ii] = losses[ii] / (2 * pi * speed[ii] / 60) * self.machine.rotor.L1

    p = self.machine.get_pole_pair_number()
    interp_fun = interp.interp1d(speed / 60 * p, t_fric, kind="cubic")

    return lambda Id, Iq, freq: interp_fun(abs(freq))
