import numpy as np


def interpolate(self, symbol, nd=21, nq=21, q_symmetry=0, machine=None):
    """
    Interpolation on datakeepers.
    """
    # get the interpolator
    interpolator, x, y = self.get_interp(symbol, q_symmetry=q_symmetry, machine=machine)

    # go linearly in the grid
    xline = np.linspace(min(x), max(x), nd)
    yline = np.linspace(min(y), max(y), nq)

    # construct 2d grid from these
    xi, yi = np.meshgrid(xline, yline)

    # interpolate z data; same shape as xgrid and ygrid
    zi = interpolator(xi, yi)

    return xi, yi, zi
