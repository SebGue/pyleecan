def get_torque_norm(self, machine):
    """Get the torque normalization factor."""
    # Fill factor and current density is assumed to be the same as in the reference

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

    xi = machine.stator.winding.comp_winding_factor()
    L1 = machine.stator.L1

    torque_norm = xi / xi_ref * L1 / L1_ref

    return torque_norm[0]
