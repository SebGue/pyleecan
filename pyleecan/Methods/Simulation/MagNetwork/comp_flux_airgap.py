# -*- coding: utf-8 -*-
def comp_flux_airgap(self, output, axes_dict, Is_val=None, Ir_val=None):
    """Build and solve FEMM model to calculate and store magnetic quantities

    Parameters
    ----------
    self : MagFEMM
        a MagFEMM object
    output : Output
        an Output object
    axes_dict: {Data}
        Dict of axes used for magnetic calculation

    Returns
    -------
    out_dict: dict
        Dict containing the following quantities:
            Br : ndarray
                Airgap radial flux density (Nt,Na) [T]
            Bt : ndarray
                Airgap tangential flux density (Nt,Na) [T]
            Tem : ndarray
                Electromagnetic torque over time (Nt,) [Nm]
            Phi_wind_stator : ndarray
                Stator winding flux (qs,Nt) [Wb]
            Phi_wind : dict
                Dict of winding fluxlinkage with respect to Machine.get_lam_list_label (qs,Nt) [Wb]
            meshsolution: MeshSolution
                MeshSolution object containing magnetic quantities B, H, mu for each time step
    """

    Bx, By, Bx_airgap, By_airgap = self.run_radial(axes_dict, Is_val)

    out_dict = {
        "B_{rad}": Bx_airgap,
        "B_{circ}": By_airgap,
    }

    return out_dict
