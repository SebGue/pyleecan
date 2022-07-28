# -*- coding: utf-8 -*-
def comp_flux_airgap(
    self, output, axes_dict, Is_val=None, Ir_val=None, type_coord_sys=2
):
    """Build and solve the MagNetwork model to calculate and store magnetic quantities

    Parameters
    ----------
    self : MagNetwork
        a MagNetwork object
    output : Output
        an Output object
    axes_dict: {Data}
        Dict of axes used for magnetic calculation
    type_coord_sys : integer (Default = 2)
        Type of the coordinate system : 1 for cartesian, 2 for Radial

    Returns
    -------
    out_dict: dict
        Dict containing the following quantities:
            Br : ndarray
                Airgap radial flux density (Nt,Na) [T]
            Bt : ndarray
                Airgap tangential flux density (Nt,Na) [T]
    """

    if type_coord_sys == 1:
        Bx, By, Bx_airgap, By_airgap = self.run_cartersian(axes_dict, Is_val)

    if type_coord_sys == 2:
        Bx, By, Bx_airgap, By_airgap = self.run_radial(axes_dict, Is_val)

    out_dict = {
        "B_{rad}": Bx_airgap,
        "B_{circ}": By_airgap,
    }

    return out_dict
