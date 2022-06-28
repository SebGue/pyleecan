# Method to compute the air gap magnetic flux density
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

# def comp_flux_airgap(self, output, axes_dict, Is_val=None, Ir_val=None):
def comp_flux_airgap_local(
    self, r, theta, F, list_elem, list_coord, la, R_ext_lam, R_int_lam
):

    # Get the airgap circular pattern
    # average_pattern = self.Machine.rotor.Rext + self.Machine.comp_width_airgap_mec() / 2

    # Compute B
    Bx, By = self.compute_B_radial(F, list_elem, list_coord, la)

    # Reshape B in 2d array
    Bx = Bx.reshape((r.size - 1, theta.size - 1))
    By = By.reshape((r.size - 1, theta.size - 1))

    # select the index of the line
    i = (int)((R_ext_lam - R_int_lam) * 1000)
    # i = 41

    # Looking for the position in the airgap
    position = r.min() + i * (r.max() - r.min()) / (r.size - 1)

    print("The horizental line in the airgap is: ", position, i)

    # Plot the line
    plt.plot(Bx[i, :], label="x conponent")
    # plt.plot(By[i, :], label="y component")

    plt.legend()
    plt.show()

    return
