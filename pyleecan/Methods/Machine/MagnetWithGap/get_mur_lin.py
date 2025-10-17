# -*- coding: utf-8 -*-


def get_mur_lin(self):
    """Get relative magnetic permeability

    Parameters
    ----------
    self : MatMagnetics
        a MatMagnetics object

    Returns
    -------
    mur_lin: float
        relative magnetic permeability

    """
    mur = self.mat_type.mur_lin

    return mur
