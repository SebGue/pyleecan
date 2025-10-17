# -*- coding: utf-8 -*-


def get_Hc(self, T_op=None, T_ref=20):
    """Get magnetic excitation coercitivity

    Parameters
    ----------
    self : MatMagnetics
        a MatMagnetics object
    T_op: float
        Material operational temperature [degC]
    T_ref: float
        Material reference temperature [degC]

    Returns
    -------
    Hc: float
        Magnetic excitation coercitivity

    """
    Hc = self.mat_type.mag.get_Hc(T_op=T_op, T_ref=T_ref)

    return Hc
