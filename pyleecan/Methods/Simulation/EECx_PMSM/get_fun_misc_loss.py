# -*- coding: utf-8 -*-
from numpy import linspace, zeros, pi, array, meshgrid
import scipy.interpolate as interp

MTH_NAME = "get_fun_core_loss()"  # TODO ?? move to __init__.py ??
CLS_NAME = "EECx_PMSM"


def get_fun_misc_loss(self, max_speed):
    """Method to get the additional losses function."""
    if self.misc_loss is None or not self.misc_loss:
        return lambda Id, Iq, freq: 0

    # some settings
    ND = 21  # number of d-current interpolates
    NQ = 21  # number of q-current interpolates
    NN = 21  # number of speed interpolates
    EXP = 2  # exponent to 'compress' loss data range, i.e. interpol on loss^(1/EXP)

    # ref. speed vector for loss calculation
    speed_ref = linspace(0, max_speed, NN).tolist()

    # ref. current for loss calculation
    sd = array(self.table["Sd"].result)
    sq = array(self.table["Sq"].result)
    kC = self.table.simu.get_current_norm(self.machine)

    Id = sd * kC
    Iq = sq * kC

    Idq = array([Id, Iq]).T

    # get the losses for all data points in the table
    losses = []

    for out in self.table.output_list:
        # store original value for speed and set desired speed
        for loss_mdl in self.misc_loss.values():
            loss = zeros(NN)
            loss_mdl = loss_mdl.copy()
            loss_mdl.N0 = speed_ref
            loss_data = loss_mdl.comp_loss(output=out)[0]
            if loss_data is None:
                data = zeros(NN)
                self.get_logger.warning()
            else:
                # get data axes to compute means and sums over the axes
                axes = loss_data.get_axes()
                axes_list = [
                    ax.name + "=sum"
                    for ax in axes
                    if ax.name not in ["time", "freqs", "speed"]
                ]

                data_dict = loss_data.get_along("time=mean", *axes_list, "speed")
                data = data_dict[loss_data.symbol]

            loss += data

        losses.append(loss.tolist())

    losses = array(losses) * self.machine.stator.L1

    # interpolate in current plane to have better resolution for linear 3D interpol.
    Idi = linspace(Id.min(), Id.max(), ND)
    Iqi = linspace(Iq.min(), Iq.max(), NQ)

    x, y = meshgrid(Idi, Iqi, indexing="ij")
    Idqi = array([x.T, y.T]).T

    zidata = zeros([ND, NQ, NN])
    for ii in range(NN):
        zdata = losses[:, ii] ** (1 / EXP)
        zidata[:, :, ii] = interp.griddata(Idq, zdata, Idqi, method="cubic")

    # interpolate all data
    p = self.machine.get_pole_pair_number()
    freq = array(speed_ref) / 60 * p
    x, y, z = meshgrid(Idi, Iqi, freq, indexing="ij")
    ipoints = array([X.flatten(order="F") for X in [x, y, z]]).T
    interp_fun = interp.LinearNDInterpolator(ipoints, zidata.flatten(order="F"))

    return lambda Id, Iq, freq: interp_fun(Id, Iq, abs(freq)) ** EXP
