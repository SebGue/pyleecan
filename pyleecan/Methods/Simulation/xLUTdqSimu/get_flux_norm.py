def get_flux_norm(self, machine):
    """Get the fluxlinkage normalization factor."""
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
    xi_ref = ref_machine.stator.winding.comp_winding_factor()
    L1_ref = ref_machine.stator.L1

    Ntspc = machine.stator.winding.comp_Ntsp()
    xi = machine.stator.winding.comp_winding_factor()
    L1 = machine.stator.L1

    flux_norm = xi / xi_ref * L1 / L1_ref * Ntspc

    return flux_norm[0]
