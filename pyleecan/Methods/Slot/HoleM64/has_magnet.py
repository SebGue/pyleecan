def has_magnet(self):
    """Return if the Hole has magnets

    Parameters
    ----------
    self : HoleM64
        A HoleM64 object

    Returns
    -------
    has_magnet : bool
        True if the magnets are not None
    """

    return self.magnet_0 is not None
