def get_torque(self, machine=None):
    """Get the torque of a (compatible) machine. If no machine is specified,
    normalized quantities will be returned.
    Paramter
    --------
    machine : a Machine object

    Return
    ------
    X, Y : narray of the d- and q-current or current density respectively
    Psid, Psiq : ndarray of the fluxlinkage
    """

    # get the data and interpolate
    Iq = self["Sq"].result

    if min(Iq) < 0:
        q_sym = 0
    else:
        q_sym = -1

    Idi, Iqi, Ti = self.interpolate(
        symbol="Tem_av", nd=101, nq=101, q_symmetry=q_sym, machine=machine
    )

    return Idi, Iqi, Ti
