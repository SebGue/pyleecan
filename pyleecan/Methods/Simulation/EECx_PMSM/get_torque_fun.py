# -*- coding: utf-8 -*-


def get_torque_fun(self):
    """Method to get a torque function T = f(Id, Iq, freq)."""
    machine = self.machine
    Ti = self.table.get_interp(symbol="Tem_av", q_symmetry=-1, machine=machine)[0]

    return lambda Id, Iq, freq: Ti(Id, Iq)
