# -*- coding: utf-8 -*-
from numpy import array, linspace, meshgrid, zeros
import scipy.interpolate as interp

MTH_NAME = "get_fun_core_loss()"
CLS_NAME = "EECx_PMSM"


def get_fun_core_loss(self, max_speed):
    """Method to get the core losses function. Therefore the given core loss models
    will be applied to the machine characteristic table to get the losses. The
    calculated losses will then be used to create a loss interpolation function.
    """
    # some settings
    ND = 21  # number of d-current interpolates
    NQ = 21  # number of q-current interpolates
    NN = 21  # number of speed interpolates
    EXP = 2  # exponent to 'compress' loss data range, i.e. interpol on loss^(1/EXP)

    # ref. speed vector for loss calculation
    speed_ref = linspace(0, max_speed, NN).tolist()

    if self.core_loss is None:
        return lambda Id, Iq, freq: 0

    keys = list(self.core_loss.keys())
    if not keys:
        return lambda Id, Iq, freq: 0

    # iterate the table data to get the losses of each pair of currents
    losses = []
    for out in self.table.output_list:
        loss = 0
        for key in keys:
            # get loss model, set speeds and calculate losses
            loss_model = self.core_loss[key].copy()
            loss_model.N0 = speed_ref
            loss_data = loss_model.comp_loss(out, key.split(".")[0])[0]

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

                # some checks to be sure that data are okay
                if data_dict["axes_dict_other"]:
                    self.get_logger().warning(
                        f"{CLS_NAME}.{MTH_NAME}: Other axes detected."
                    )
                if not data_dict["axes_list"] or "speed" not in data_dict.keys():
                    self.get_logger().warning(
                        f"{CLS_NAME}.{MTH_NAME}: No speed axis found."
                    )

                if len(data.shape) > 2 or (
                    len(data.shape) == 2 and 1 not in data.shape
                ):
                    self.get_logger().warning(
                        f"{CLS_NAME}.{MTH_NAME}: "
                        + "Loss data must be a vector but is of shape {data.shape}."
                    )

                if (
                    speed.tolist() != speed_ref
                ):  # TODO better handling, i.e. extend single values to speed_ref
                    self.get_logger().warning(
                        f"{CLS_NAME}.{MTH_NAME}: Data speed axis is wrong. "
                        + f"It is {speed} but should be {speed_ref}."
                    )

                loss += data * self.machine.get_lam_by_label(key.split(".")[0]).L1
        losses.append(loss.tolist())

    losses = array(losses)

    # === interpolate the loss data ===
    sd = array(self.table["Sd"].result)
    sq = array(self.table["Sq"].result)
    kC = self.table.simu.get_current_norm(self.machine)

    Id = sd * kC
    Iq = sq * kC

    Idi = linspace(Id.min(), Id.max(), ND)
    Iqi = linspace(Iq.min(), Iq.max(), NQ)

    Idq = array([Id, Iq]).T
    x, y = meshgrid(Idi, Iqi, indexing="ij")
    Idqi = array([x.T, y.T]).T

    # interpolate in current plane to have better resolution for linear 3D interpol.
    zidata = zeros([Idi.size, Iqi.size, losses.shape[1]])
    for ii in range(losses.shape[1]):
        zdata = losses[:, ii] ** (1 / EXP)
        zidata[:, :, ii] = interp.griddata(Idq, zdata, Idqi, method="cubic")
        # interpolator = interp.CloughTocher2DInterpolator(Idq, zdata)
        # zidata = interpolator(Idi, Iqi)

    # interpolate all data
    p = self.machine.get_pole_pair_number()
    freq = speed / 60 * p
    x, y, z = meshgrid(Idi, Iqi, freq, indexing="ij")
    ipoints = array([X.flatten(order="F") for X in [x, y, z]]).T
    interp_fun = interp.LinearNDInterpolator(ipoints, zidata.flatten(order="F"))

    return lambda Id, Iq, freq: interp_fun(Id, Iq, abs(freq)) ** EXP
