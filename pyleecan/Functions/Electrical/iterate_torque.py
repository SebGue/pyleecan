# -*- coding: utf-8 -*-

from scipy.optimize import minimize
from numpy import sqrt
from pprint import pprint

# TODO some variable to be moved to Electrical class
is_update_parameters = False  # update parameters to estimated OP
max_current = 1000  # [A] max. current magnitude
max_voltage = 1000  # [V] max. voltage magnitude


def _comp_voltage(I, elec, output):
    #
    param = elec.eec.parameters

    # delete voltages to force voltage computation
    if "Ud" in elec.eec.parameters:
        del elec.eec.parameters["Ud"]
    if "Uq" in elec.eec.parameters:
        del elec.eec.parameters["Uq"]

    # set currents
    param["Id"] = I[0]
    param["Iq"] = I[1]
    output.elec.Id_ref = I[0]
    output.elec.Iq_ref = I[1]

    # adapt flux
    if "Ld" in param and param["Ld"]:
        param["Phid"] = param["phi"] + I[0] * param["Ld"]

    if "Lq" in param and param["Lq"]:
        param["Phiq"] = I[1] * param["Lq"]

    # Solve the electrical equivalent circuit
    elec.eec.solve_EEC(output)

    # return the abs. voltage
    U = sqrt(output.elec.Ud_ref ** 2 + output.elec.Uq_ref ** 2)
    # print(f'U = {U}, Id = {I[0]}, Iq = {I[1]}')
    return U


def _comp_torque(I, elec, output):
    # store reference torque, since it will be overridden by comp_torque()
    torque_ref = output.elec.Tem_av_ref

    # update currents and voltages
    _comp_voltage(I, elec, output)

    # compute losses, power and torque
    elec.eec.comp_joule_losses(output)
    elec.comp_power(output)
    elec.comp_torque(output)

    torque = output.elec.Tem_av_ref
    output.elec.Tem_av_ref = torque_ref

    Id = output.elec.Id_ref
    Iq = output.elec.Iq_ref
    Ud = output.elec.Ud_ref
    Uq = output.elec.Uq_ref

    # print(f'T = {torque}, Ud = {Ud}, Uq = {Uq}, Id = {Id}, Iq = {Iq}')

    return torque


def iterate_torque(elec, output):
    # Compute first set of parameters of the electrical equivalent circuit
    elec.eec.comp_parameters(output)

    pprint(elec.eec.parameters)
    print(output.elec.Id_ref, output.elec.Iq_ref)

    # Solve the electrical equivalent circuit
    # elec.eec.solve_EEC(output)

    # get reference torque
    torque_ref = output.elec.Tem_av_ref

    # tolerated difference in torque computation
    torque_err = torque_ref / 1e4 + 1e-4  # 0.1 % + 0.1 mNm

    # objective
    objectiv = lambda x: sqrt(x[0] ** 2 + x[1] ** 2)  # minimize the current

    # boundaries
    bnds = None  # ((-max_current, max_current), (-max_current, max_current))

    # constraints
    max_current_constrain = lambda x: -sqrt(x[0] ** 2 + x[1] ** 2) + max_current
    max_voltage_constrain = lambda x: -_comp_voltage(x, elec, output) + max_voltage
    torque_constrain = (
        lambda x: -((_comp_torque(x, elec, output) - torque_ref) ** 2) + torque_err ** 2
    )

    cons = (
        # {"type": "ineq", "fun": max_current_constrain},
        # {"type": "ineq", "fun": max_voltage_constrain},
        {"type": "ineq", "fun": torque_constrain},
    )

    x_start = [10, 10]

    # res = minimize(objectiv, x_start, method='trust-constr', constraints = cons, bounds = bnds) #slower than  SLSQP
    res = minimize(objectiv, x_start, method="SLSQP", constraints=cons, bounds=bnds)

    if res.success == False:
        print("False")
        output.elec.Id_ref = None
        output.elec.Iq_ref = None
        output.elec.Is = None
        output.elec.Ud_ref = None
        output.elec.Uq_ref = None
        output.elec.Us = None
    else:
        # final solve
        # _comp_voltage(res.x, elec, output)
        pass
