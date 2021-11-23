from numpy import array
from ....Functions.Electrical.dqh_transformation import n2dq
from ....Classes.DataKeeper import DataKeeper


def _set_fluxlinkage(self):
    """Compute the normalized dq fluxlinkage and set it as a datakeeper"""
    p = self.simu.machine.get_pole_pair_number()
    qs = self.simu.machine.stator.winding.qs
    Ntspc = self.simu.machine.stator.winding.comp_Ntsp()

    theta = self.get_angle_rotor() - self.get_angle_offset_initial()
    Psi = []
    for psi in self["Psi"].result:
        # TODO use harmonics instead
        # val = psi.get_along("time", "phase")[psi.symbol]
        val = psi
        Psi_dq = n2dq(val, theta * p, n=qs, is_dq_rms=True)
        Psi.append(Psi_dq.mean(axis=0))

    Psi = Psi / Ntspc

    # Construct new datakeeper and store results
    symbol = "Psi_d"
    self.xoutput_dict[symbol] = DataKeeper(
        name="norm. dFluxlinkage",
        symbol=symbol,
        unit="Vs/m/turns",
        keeper=None,
        error_keeper=None,
        result=array(Psi)[:, 0].tolist(),
    )

    symbol = "Psi_q"
    self.xoutput_dict[symbol] = DataKeeper(
        name="norm. qFluxlinkage",
        symbol=symbol,
        unit="Vs/m/turns",
        keeper=None,
        error_keeper=None,
        result=array(Psi)[:, 1].tolist(),
    )
