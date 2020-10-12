# -*- coding: utf-8 -*-

from ....Functions.Electrical.comp_fluxlinkage import comp_fluxlinkage as comp_flx
from numpy import mean


def comp_inductance(self, output):
    """Compute using FEMM the inductance

    Parameters
    ----------
    self : IndMagFEMM
        an IndMagFEMM object
    output : Output
        an Output object
    """
    # store orignal currents
    Is = output.elec.Is
    Id_ref = output.elec.Id_ref
    Iq_ref = output.elec.Iq_ref

    if Id_ref == 0 or Iq_ref == 0:
        # Set currents to 1A for the FEMM simulation to get Ld, Lq parameters
        output.elec.Is = None
        output.elec.Id_ref = 1  # TODO find resonable value for actual machine
        output.elec.Iq_ref = 1

    # compute the fluxlinkage
    fluxdq = comp_flx(self, output)

    # restore orignal currents
    output.elec.Is = Is
    output.elec.Id_ref = Id_ref
    output.elec.Iq_ref = Iq_ref

    return (mean(fluxdq[0]), mean(fluxdq[1]))
