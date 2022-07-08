# Method to compute the air gap magnetic flux density localy
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

# def comp_flux_airgap(self, output, axes_dict, Is_val=None, Ir_val=None):
def comp_flux_airgap_local(self, r, theta, F, list_elem, list_coord, la, Rgap):

    Machine = self.parent.machine
    # Compute B
    Bx, By = self.compute_B_radial(F, list_elem, list_coord, la)

    # Reshape B in 2d array
    Bx = Bx.reshape((r.size - 1, theta.size - 1))
    By = By.reshape((r.size - 1, theta.size - 1))

    # select the index of the line
    i = (int)((Machine.stator.Rint - Machine.rotor.Rint) * 1000)
    # i = (int)(Rgap * 1000)

    # Looking for the position in the airgap
    position = r.min() + i * (r.max() - r.min()) / (r.size - 1)

    print("The horizental line in the airgap is: ", position, i)

    # Plot the line
    # plt.plot(Bx[i, :], label="x conponent")
    # plt.plot(By[i, :], label="y component")

    # plt.legend()
    # plt.show()

    return Bx[i, :], By[i, :]
