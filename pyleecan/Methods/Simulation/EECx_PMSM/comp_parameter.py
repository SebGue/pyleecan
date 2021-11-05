# -*- coding: utf-8 -*-


def comp_parameter(self, max_speed):
    # === get the losses functions ===
    comp_loss_core = self.get_fun_core_loss(max_speed)
    comp_loss_misc = self.get_fun_misc_loss(max_speed)
    comp_loss_mech = self.get_fun_mech_loss(max_speed)
    comp_Rac = self.winding.get_resistance_fun(self.machine, self.tW)

    # === get fluxlinkage and torque functions ====
    comp_torque = self.get_torque_fun()
    comp_fluxlinkage = self.get_fluxlinkage_fun()

    # === setup model parameters ===
    pmsm_param = dict()
    pmsm_param["TorqueFunc"] = comp_torque
    pmsm_param["FluxlinkageFunc"] = comp_fluxlinkage
    pmsm_param["CoreLossFunc"] = comp_loss_core
    pmsm_param["MechLossFunc"] = comp_loss_mech
    pmsm_param["MiscLossFunc"] = comp_loss_misc
    pmsm_param["StatorResistanceFunc"] = comp_Rac

    # === setup parameter space bounds ===
    Psidi, Id, Iq = self.table.get_interp(
        symbol="Psi_d", machine=self.machine, q_symmetry=0
    )

    pmsm_param["Idmin"] = min(Id)
    pmsm_param["Idmax"] = max(Id)
    pmsm_param["Iqmin"] = min(Iq)
    pmsm_param["Iqmax"] = max(Iq)

    return pmsm_param
