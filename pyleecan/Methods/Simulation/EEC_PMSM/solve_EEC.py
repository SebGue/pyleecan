# -*- coding: utf-8 -*-

from numpy import array, pi
from scipy.linalg import solve


def solve_EEC(self, output):
    """Compute the parameters dict for the equivalent electrical circuit
    cf "Advanced Electrical Drives, analysis, modeling, control"
    Rik de doncker, Duco W.J. Pulle, Andre Veltman, Springer edition

                 <---                               --->
     -----R-----wsLqIq----              -----R-----wsLdId----
    |                     |            |                     |
    |                     |            |                    BEMF
    |                     |            |                     |
     ---------Id----------              ---------Iq----------

             --->                               --->
              Ud                                 Uq

    Parameters
    ----------
    self : EEC_PMSM
        an EEC_PMSM object
    output : Output
        an Output object
    """
    # readability
    felec = output.elec.felec
    ws = 2 * pi * felec
    R = self.parameters["R20"]
    Ld = self.parameters["Ld"]
    Lq = self.parameters["Lq"]

    Phi = self.parameters["Phi"]

    # Prepare linear system

    # Solve system
    if "Ud" in self.parameters:
        Ud = self.parameters["Ud"]
        Uq = self.parameters["Uq"]

        XR = array(
            [
                [R, -ws * Lq],
                [ws * Ld, R],
            ]
        )
        XE = array([0, ws * Phi])
        XU = array([Ud, Uq])
        XI = solve(XR, XU - XE)
        output.elec.Id_ref = XI[0]
        output.elec.Iq_ref = XI[1]
    else:
        Id = self.parameters["Id"]
        Iq = self.parameters["Iq"]
        output.elec.Ud_ref = R * Id - ws * Lq * Iq
        output.elec.Uq_ref = R * Iq + ws * Ld * Id + ws * Phi

    # Compute currents
    output.elec.Is = None
    output.elec.Is = output.elec.get_Is()

    # Compute voltage
    output.elec.Us = None
    output.elec.Us = output.elec.get_Us()
