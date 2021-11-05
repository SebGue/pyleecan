from ....Classes.CondType22 import CondType22


def get_norm_machine(self, machine=None):
    """Copy an normalize a machine to conform to the MagLUTSimu.
    Therefore the lenght is set to 1m, parallel paths are set to 1, set magnet
    temperature and the conductors are substituted with a generic conductor with
    the specified fill factor. All other winding parameters will stay the same.
    """
    # copy machine first
    machine = machine.copy()

    # set length and parallel paths
    machine.stator.winding.Npcp = 1
    machine.stator.L1 = 1
    machine.rotor.L1 = 1

    # set generic conductor
    Ncps = _comp_Ncps(machine.stator.winding)
    S_slot_wind = machine.stator.slot.comp_surface_active()

    cond_mat = machine.stator.winding.conductor.cond_mat
    ins_mat = machine.stator.winding.conductor.ins_mat

    machine.stator.winding.conductor = CondType22()
    machine.stator.winding.conductor.cond_mat = cond_mat
    machine.stator.winding.conductor.ins_mat = ins_mat
    machine.stator.winding.conductor.Sbar = S_slot_wind / Ncps * self.fill_factor

    # set magnet temperature
    magnets = machine.rotor.get_magnet_list()
    for magnet in magnets:
        magnet.set_temp_magnet(self.tM)

    return machine


# TODO move to LamSlotWind or use comp.fill_factor()
def _comp_Ncps(winding):
    # compute the number of conductors per slot
    Ncps_ = abs(winding.get_connection_mat().sum(axis=(0, 1))).sum(axis=1)
    Ncps = Ncps_.mean()

    if Ncps_.std() != 0:
        winding.get_logger().warning(
            "LamSlotWind.comp_fill_factor: "
            "Uneven number of conductors per slot. "
            + "Max. number of conductors will be returned."
        )
        Ncps = Ncps_.max()

    return Ncps
