# -*- coding: utf-8 -*-


def get_rho_elec(self):
    """Get electrical resistivity

    Parameters
    ----------
    self : MatMagnetics
        a MatMagnetics object

    Returns
    -------
    rho: float
        electrical resistivity

    """
    rho = self.mat_type.elec.rho

    return rho
