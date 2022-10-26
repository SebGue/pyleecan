# -*- coding: utf-8 -*->


def build_geometry(self, alpha=0, delta=0, is_simplified=False):
    """Abstract class method

    Parameters
    ----------
    self : Hole
        A Hole object
    alpha : float
        Angle to rotate the slot (Default value = 0) [rad]
    delta : complex
        Complex to translate the slot (Default value = 0)
    is_simplified : bool
       True to avoid line superposition

    Returns
    -------
    surf_list: list
        List of SurfLine needed to draw the HoleM50

    """

    return []
