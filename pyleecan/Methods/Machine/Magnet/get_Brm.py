# -*- coding: utf-8 -*-


def get_Brm(self, T_op=None, T_ref=20):
    """Get magnetic remanent flux density

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
    Brm: float
        Magnetic remanent flux density

    """
    Brm = self.mat_type.mag.get_Brm(T_op=T_op, T_ref=T_ref)

    return Brm
