def get_current_norm(self, machine):
    """Get the current normalization factor (including current density)."""
    # check if machine is feasable
    if not self.is_norm_machine(machine):
        self.get_logger().warning("MagLUTSimu: Given machine is not feasable.")
        return None

    if machine.rotor.L1 != machine.stator.L1:
        self.get_logger().warning(
            "MagLUTSimu: Given machine has different stator and rotor lenght."
        )
        return None

    # get the normalized machine
    ref_machine = self.get_norm_machine(self.machine)
    Ntspc_ref = ref_machine.stator.winding.comp_Ntsp()
    Npcp_ref = ref_machine.stator.winding.Npcp
    S_cond_ref = ref_machine.stator.winding.conductor.comp_surface_active()

    Ntspc = machine.stator.winding.comp_Ntsp()
    Npcp = machine.stator.winding.Npcp

    # equations
    # I_ref = Dens_ref * S_cond_ref * Npcp_ref
    # I = Dens * S_cond * Npcp
    # I = I_ref * Ntspc_ref / Ntspc
    #
    # I = Dens_ref * S_cond_ref * Npcp_ref * Ntspc_ref / Ntspc

    current_norm = S_cond_ref * Npcp_ref * Ntspc_ref / Ntspc

    return current_norm
