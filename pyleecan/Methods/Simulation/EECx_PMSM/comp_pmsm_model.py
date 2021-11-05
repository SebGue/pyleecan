# -*- coding: utf-8 -*-
from numpy import array, pi

# =============================================================================
# machine modell with current and frequency input
def comp_pmsm_model(self, Id, Iq, speed, parameter_dict):
    """Compute the PMSM model quantities."""
    # TODO regard fill factor and skew factor

    # get input parameter from dict
    p = self.machine.stator.get_pole_pair_number()
    qs = self.machine.stator.winding.qs

    fTorque = parameter_dict["TorqueFunc"]
    fPsi = parameter_dict["FluxlinkageFunc"]
    fPvu = parameter_dict["CoreLossFunc"]
    fPadd = parameter_dict["MiscLossFunc"]
    fTfriction = parameter_dict["MechLossFunc"]
    fR = parameter_dict["StatorResistanceFunc"]

    # ---
    freq = speed / 60 * p

    # calculate torque and fluxlinkage
    Psi_d, Psi_q = fPsi(Id, Iq, freq)
    Tem = 3 * p * (Psi_d * Iq - Psi_q * Id)  # fTorque(Id, Iq, freq)
    Tfric = fTfriction(Id, Iq, freq)

    # shaft torque
    if speed > 0:
        Tshaft = Tem - abs(Tfric)
    elif speed < 0:
        Tshaft = Tem + abs(Tfric)
    elif speed == 0:
        Tshaft = Tem  # TODO check this case

    # calculate iron and misc. losses
    Pvu = fPvu(Id, Iq, freq)
    Padd = fPadd(Id, Iq, freq)

    # calculate induced voltage
    omega = 2 * pi * freq

    Uxd = -omega * Psi_q
    Uxq = omega * Psi_d

    Ux_sq = Uxd ** 2 + Uxq ** 2

    # calculate equivalent iron loss resistance and current
    if Pvu < 1e-3:
        Ife_d = 0
        Ife_q = 0
    else:
        Rfe = qs * Ux_sq / Pvu
        Ife_d = Uxd / Rfe
        Ife_q = Uxq / Rfe

    # corrected current and resistive voltage drop
    Id_ = Id + Ife_d
    Iq_ = Iq + Ife_q
    Rac = fR(Id_, Iq_, freq)

    Urd = Id_ * Rac
    Urq = Iq_ * Rac

    # losses
    Pvcu = qs * (Urd * Id_ + Urq * Iq_)
    Pvmech = 2 * pi * speed / 60 * Tfric

    # armature voltage
    Ud = Urd + Uxd
    Uq = Urq + Uxq

    return Ud, Uq, Id_, Iq_, Tshaft, Pvu, Pvcu, Pvmech, Padd, Rac
