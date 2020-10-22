# -*- coding: utf-8 -*-


def comp_power(self, output):
    """Compute the electrical average power

    Parameters
    ----------
    self : Electrical
        an Electrical object
    output : Output
        an Output object
    """

    qs = output.simu.machine.stator.winding.qs
    Id = output.elec.Id_ref
    Iq = output.elec.Iq_ref
    Ud = output.elec.Ud_ref
    Uq = output.elec.Uq_ref

    if Ud is None or Uq is None or Id is None or Iq is None:
        Pem_av_ref = None
    else:
        Pem_av_ref = qs * (Ud * Id + Uq * Iq) / 2

    output.elec.Pem_av_ref = Pem_av_ref
