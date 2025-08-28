def get_magnet_dict(self):
    """Return a dictionary with all the magnets of the Hole

    Parameters
    ----------
    self : Hole
        A Hole object

    Returns
    -------
    magnet_dict : {Magnet}
        Dictionary of magnet (key = magnet_X, value= Magnet or None)
    """

    return dict(
        magnet_0=self.magnet_0,
    )
