# -*- coding: utf-8 -*-

from scipy.optimize import minimize
from numpy import sqrt
from pprint import pprint

# TODO some variable to be moved to Electrical class
is_update_parameters = False  # update parameters to estimated OP
max_current = 5000  # [A] max. current magnitude
max_voltage = 563.38  # [V] max. voltage magnitude
TOL = 1e-3  # [] aceptable tolerance of current difference


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
    output.elec.Is = None
    output.elec.Id_ref = I[0]
    output.elec.Iq_ref = I[1]

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

    # pprint(elec.eec.parameters)
    # print(output.elec.Id_ref, output.elec.Iq_ref)

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
    max_current_constrain = lambda x: -sqrt(x[0] ** 2 + x[1] ** 2) + elec.Imax
    max_voltage_constrain = lambda x: -_comp_voltage(x, elec, output) + elec.Umax
    torque_constrain = (
        lambda x: -((_comp_torque(x, elec, output) - torque_ref) ** 2) + torque_err ** 2
    )

    cons = [
        {"type": "ineq", "fun": torque_constrain},
        {"type": "ineq", "fun": max_current_constrain},
        {"type": "ineq", "fun": max_voltage_constrain},
    ]

    if elec.Imax is None:
        cons.pop(1)
    if elec.Umax is None:
        cons.pop(-1)

    cons = tuple(cons)

    x_start = [10, 10]

    # res = minimize(objectiv, x_start, method='trust-constr', constraints = cons, bounds = bnds) #slower than  SLSQP
    res = minimize(objectiv, x_start, method="SLSQP", constraints=cons, bounds=bnds)

    if res.success == True:
        # final solve needed ???
        # _comp_voltage(res.x, elec, output)
        if is_update_parameters:
            pass
            """ Does not work ATM
            # store actual I_ref values
            Id_ref = output.elec.Id_ref
            Iq_ref = output.elec.Iq_ref
            I = abs(Id_ref + 1j * Iq_ref)

            DO_ITER = True
            i = 0
            while DO_ITER:  # TODO test parameter diff. instead ??? maybe not good
                i += 1
                print("Parameter Iteration " + str(i))
                # recompute parameters
                # TODO improve by comp. diff. ind. around OP -> iterate over diff ind first -> have Ld/Lq adaption function
                elec.eec.parameters["Id"] = None  # to force recomp.
                elec.eec.parameters["Iq"] = None
                # print(output.elec.Id_ref, output.elec.Iq_ref)
                elec.eec.comp_parameters(output)
                pprint(elec.eec.parameters)
                print(output.elec.Id_ref, output.elec.Iq_ref)

                # recompute electrical quantities
                x_start = res.x
                res = minimize(
                    objectiv, x_start, method="SLSQP", constraints=cons, bounds=bnds
                )

                print(output.elec.Id_ref, output.elec.Iq_ref)

                # print(res.x)
                if res.success == False:
                    DO_ITER = False
                else:
                    Id_ref_ = output.elec.Id_ref
                    Iq_ref_ = output.elec.Iq_ref
                    cond_1 = abs(Id_ref - Id_ref_) <= (TOL * abs(Id_ref) + TOL * I)
                    cond_2 = abs(Iq_ref - Iq_ref_) <= (TOL * abs(Iq_ref) + TOL * I)
                    # print(cond_1, cond_2)
                    # print(Id_ref, Id_ref_, Iq_ref, Iq_ref_)
                    if cond_1 and cond_2:
                        DO_ITER = False
                    else:
                        # adapt for (more) stable iteration
                        k = 0.1
                        output.elec.Id_ref = k * Id_ref_ + (1 - k) * Id_ref
                        output.elec.Iq_ref = k * Iq_ref_ + (1 - k) * Iq_ref

                        Id_ref = Id_ref_
                        Iq_ref = Iq_ref_
                        I = abs(Id_ref + 1j * Iq_ref)
            """

        output.elec.Is = None
        output.elec.Is = output.elec.get_Is()

    if res.success == False:
        print("False")
        output.elec.Id_ref = None
        output.elec.Iq_ref = None
        output.elec.Is = None
        output.elec.Ud_ref = None
        output.elec.Uq_ref = None
        output.elec.Us = None

