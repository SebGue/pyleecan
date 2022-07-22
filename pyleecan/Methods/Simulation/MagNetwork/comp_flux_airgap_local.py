# Method to compute the air gap magnetic flux density localy
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

# def comp_flux_airgap(self, output, axes_dict, Is_val=None, Ir_val=None):
def comp_flux_airgap_local(self, r, theta, F, list_elem, list_coord, la, Rgap):

    # Getting the machine object
    Machine = self.parent.machine

    # Computing the radial and tangential flux density
    Bx, By = self.compute_B_radial(F, list_elem, list_coord, la)

    # Reshaping B in a 2D array
    Bx = Bx.reshape((r.size - 1, theta.size - 1))
    By = By.reshape((r.size - 1, theta.size - 1))

    # Looking for the position in the center of the mecanical airgap
    position = Rgap

    # Getting the index of the line i from the position in the center of the mecanical airgap
    i = (int)((position - r.min()) * ((r.size - 1) / (r.max() - r.min())))

    # print("The horizental line in the airgap is: ", position, i)

    return Bx[i, :], By[i, :]
