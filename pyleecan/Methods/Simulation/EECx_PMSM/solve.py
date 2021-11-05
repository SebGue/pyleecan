# -*- coding: utf-8 -*-
from numpy import full, nan, ndindex, abs, sqrt, pi
from numpy.core.numeric import zeros_like
from scipy.optimize import minimize


def solve(self, ref_data, param):
    """Method to solve the EEC of a PMSM on one set of operation points."""
    # === get some values for readability ===
    speed = ref_data["speed"]
    torque = ref_data["torque"]
    mask = ref_data["mask"]

    param["Umax"] = ref_data["voltage_max"]
    if param["Umax"] is not None:  # model uses phase values
        param["Umax"] = param["Umax"] / sqrt(3)
    param["Imax"] = ref_data["current_max"]

    # === setup results arrays ===
    siz = speed.shape
    Ud, Uq = full(siz, nan), full(siz, nan)
    Id, Iq = full(siz, nan), full(siz, nan)
    Id_int, Iq_int = full(siz, nan), full(siz, nan)
    Torque = full(siz, nan)
    Pvu, Pvcu, Pvmech = full(siz, nan), full(siz, nan), full(siz, nan)
    Padd = full(siz, nan)
    Rac = full(siz, nan)

    torque0 = torque[0] if torque.ndim == 1 else torque[0, 0]
    if torque0 == 0:
        Idq0 = None
    elif torque0 > 0:
        Idq0 = [param["Idmin"] / 2, param["Iqmax"] / 2]
    else:
        Idq0 = [param["Idmin"] / 2, -param["Iqmax"] / 2]

    # === solve ===
    for idx in ndindex(speed.shape):
        if mask is None or mask[idx] == 0:
            # tolerated difference in torque computation
            TRQ_ERR = abs(torque[idx] / 1e4 + 1e-4)

            # solve
            res = _minimize_SLSQP(self, Idq0, speed[idx], torque[idx], TRQ_ERR, param)
            if res.success == True:
                Idq0 = res.x  # start value for next OP
                out_args = self.comp_pmsm_model(*res.x, speed[idx], param)

                Ud[idx], Uq[idx] = out_args[0:2]
                Id[idx], Iq[idx] = out_args[2:4]
                Torque[idx] = out_args[4]
                Pvu[idx], Pvcu[idx], Pvmech[idx], Padd[idx] = out_args[5:9]
                Rac[idx] = out_args[9]
                Id_int[idx], Iq_int[idx] = res.x[0], res.x[1]
            else:
                pass  # breakpoint

    U = (Ud ** 2 + Uq ** 2) ** (1 / 2) * sqrt(3)
    I = (Id ** 2 + Iq ** 2) ** (1 / 2)
    S = I * U * sqrt(3)
    Pel = (Ud * Id + Uq * Iq) * 3

    Pv = Pvu + Pvcu + Pvmech + Padd

    power_factor = zeros_like(S)
    ii = S != 0
    power_factor[ii] = Pel[ii] / S[ii]
    power_factor.clip(-1, 1)

    # mech. power from FEA airgap torque; may differ from ... psi * i
    Pmech = 2 * pi * speed / 60 * Torque

    Efficiency = zeros_like(Pmech)
    ii = Pmech > 0
    Efficiency[ii] = Pmech[ii] / (Pmech[ii] + Pv[ii])
    ii = Pmech < 0
    Efficiency[ii] = Pel[ii] / (Pel[ii] + Pv[ii])

    Acond = self.machine.stator.winding.conductor.comp_surface_active()
    Npcp = self.machine.stator.winding.Npcp
    J = I / (Acond * Npcp)

    res_data = dict(
        U=U,
        I=I,
        Ud=Ud * sqrt(3),
        Uq=Uq * sqrt(3),
        Id=Id,
        Iq=Iq,
        Idint=Id_int,
        Iqint=Iq_int,
        T=Torque,
        Pel=Pel,
        Pmech=Pmech,
        Pv=Pv,
        Pvu=Pvu,
        Pvcu=Pvcu,
        Pvmech=Pvmech,
        Padd=Padd,
        Rac=Rac,
        PF=power_factor,
        S=S,
        J=J,
        eta=Efficiency,
    )

    return res_data


# === helper functions ===
def comp_voltage(self, Id, Iq, speed, param):
    Ud, Uq = self.comp_pmsm_model(Id, Iq, speed, param)[0:2]
    U = (Ud ** 2 + Uq ** 2) ** (1 / 2)
    return U, Ud, Uq


def comp_current(self, Id, Iq, speed, param):
    Id_, Iq_, = self.comp_pmsm_model(
        Id, Iq, speed, param
    )[2:4]
    I = (Id_ ** 2 + Iq_ ** 2) ** (1 / 2)
    return I, Id_, Iq_


def comp_trq_diff(self, Id, Iq, speed, param, torque):
    return abs(self.comp_pmsm_model(Id, Iq, speed, param)[4] - torque)


def comp_friction(self, *args):
    return 0


def _minimize_SLSQP(self, Idq0, speed, torque, TRQ_ERR, param):
    # TODO maybe replace by cvxpy or cobyla?
    # === set bounds ===
    # bounds due to parameter space
    Id_bnds = (param["Idmin"], param["Idmax"])
    Iq_bnds = (param["Iqmin"], param["Iqmax"])
    bnds = (Id_bnds, Iq_bnds)

    # drive bounds
    Imax = param["Imax"]
    Umax = param["Umax"]

    # start value
    if Idq0 is None or None in Idq0:
        Idq0 = [-0.1, 1]

    # constraints
    EXP = 1
    torque_con = lambda x: TRQ_ERR ** EXP - (
        comp_trq_diff(self, *x, speed, param, torque) ** EXP
    )

    cons = []
    cons.append({"type": "ineq", "fun": torque_con})
    if Imax is not None:
        max_current_cons = lambda x: Imax - comp_current(self, *x, speed, param)[0]
        cons.append({"type": "ineq", "fun": max_current_cons})
    if Umax is not None:
        max_voltage_con = lambda x: Umax - comp_voltage(self, *x, speed, param)[0]
        cons.append({"type": "ineq", "fun": max_voltage_con})
    cons = tuple(cons)

    # objective: minimize the current
    objectiv = lambda x: comp_current(self, *x, speed, param)[0]

    res = minimize(objectiv, Idq0, method="SLSQP", constraints=cons, bounds=bnds)

    return res
