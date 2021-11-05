def get_fluxlinkage(self, machine=None):
    """Get the fluxlinkage of a (compatible) machine. If no machine is specified,
    normalized quantities will be returned.
    Paramter
    --------
    machine : a Machine object

    Return
    ------
    X, Y : narray of the d- and q-current or current density respectively
    Psid, Psiq : ndarray of the fluxlinkage
    """

    # set symmetry
    Iq = self["Sq"].result
    q_sym = 0 if min(Iq) < 0 else -1

    # interpolate
    Idi, Iqi, Pdi = self.interpolate(
        symbol="Psi_d", nd=101, nq=101, q_symmetry=1, machine=machine
    )

    Idi, Iqi, Pqi = self.interpolate(
        symbol="Psi_q", nd=101, nq=101, q_symmetry=q_sym, machine=machine
    )

    return Idi, Iqi, Pdi, Pqi
