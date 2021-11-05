# -*- coding: utf-8 -*-


def get_fluxlinkage_fun(self):
    """Method to get a fluxlinkage funciton Psid, Psiq = f(Id, Iq, freq)."""
    machine = self.machine
    Psidi = self.table.get_interp(symbol="Psi_d", q_symmetry=0, machine=machine)[0]
    Psiqi = self.table.get_interp(symbol="Psi_q", q_symmetry=-1, machine=machine)[0]

    return lambda Id, Iq, freq: (Psidi(Id, Iq), Psiqi(Id, Iq))
