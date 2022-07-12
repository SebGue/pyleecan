# -*- coding: utf-8 -*-


def currents(Is_val=1):

    # Machine = self.parent.machine

    # nb_phases = Machine.stator.winding.qs
    # surface_active = self.parent.machine.stator.slot.comp_surface_active()
    nb_phases = 2
    surface_active = 0.00002

    current = []
    j = []

    if Is_val is not None:
        for i in range(nb_phases):
            # j[i] = Is_val[i, :] / surface_active
            j = Is_val / surface_active
            current.append(j)
    else:
        j = None
        current.append(j)

    return current


currents(Is_val=1)
