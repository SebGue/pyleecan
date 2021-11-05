from numpy import gradient


def comp_grad_flux(self, nd=21, nq=21):
    # get the interpolated values
    xi, yi, Psi_d = self.interpolate(symbol="Psi_d", nd=nd, nq=nq, q_symmetry=1)
    xi, yi, Psi_q = self.interpolate(symbol="Psi_q", nd=nd, nq=nq, q_symmetry=-1)

    # comp. gradient
    dx = xi[0, 1] - xi[0, 0]
    dy = yi[1, 0] - yi[0, 0]

    dPsid = gradient(Psi_d, dy, dx)
    Ldd = dPsid[1]
    Ldq = dPsid[0]

    dPsiq = gradient(Psi_q, dy, dx)
    Lqd = dPsiq[1]
    Lqq = dPsiq[0]

    return xi, yi, Ldd, Ldq, Lqd, Lqq
