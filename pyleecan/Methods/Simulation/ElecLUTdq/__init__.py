# -*- coding: utf-8 -*-

from ....Classes.EEC_PMSM import EEC_PMSM

from dataclasses import dataclass
from numpy import zeros_like


@dataclass(frozen=True)
class _SIMU:
    loss: None = None


@dataclass(frozen=True)
class EEC2LUT:
    eec: EEC_PMSM
    simu: _SIMU = _SIMU()

    def interp_Phi_dqh(self, Id, Iq):
        Phi_d = Id * self.eec.Ld + self.eec.Phid_mag
        Phi_q = Iq * self.eec.Lq + self.eec.Phiq_mag
        Phi_h = zeros_like(Id)
        return Phi_d, Phi_q, Phi_h

    def interp_Tem_rip_dqh(self, Id, Iq):
        return zeros_like(Id)

    def get_Phi_dqh_mag_mean(self):
        return self.eec.Phid_mag, self.eec.Phiq_mag, 0
