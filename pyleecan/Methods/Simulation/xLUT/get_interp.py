import numpy as np
import scipy.interpolate as interp


def get_interp(self, symbol, q_symmetry=0, machine=None):
    """
    Interpolation on datakeepers.
    """
    x_symbol = "Sd"
    y_symbol = "Sq"
    z_symbol = symbol

    # check input
    symbols = [x_symbol, y_symbol, z_symbol]
    for symbol in symbols:
        if symbol not in self.keys():
            self.get_logger().warning(f"Symbol '{symbol}' not in Output.")
            return None

    # get data
    x = self[x_symbol].result.copy()
    y = self[y_symbol].result.copy()
    z = self[z_symbol].result.copy()

    # symmetry
    if q_symmetry != 0:
        x.extend([val for val in x])
        y.extend([-val for val in y])
        z.extend([np.sign(q_symmetry) * val for val in z])
        if q_symmetry == -1:
            z = [zz if yy != 0 else 0 for zz, yy in zip(z, y)]

    # return normalized quantities if no valid machine is specified
    if machine is not None:
        # check if machine is feasable, remove winding since it can differ
        if not self.simu.is_norm_machine(machine):
            self.get_logger().warning("CCOutput: Given machine is not feasable.")
            return None, None, None, None

        kC = self.simu.get_current_norm(machine)
        x = [xx * kC for xx in x]
        y = [yy * kC for yy in y]

        if "Tem_av" in z_symbol:
            kZ = self.simu.get_torque_norm(machine)
        elif "Psi" in z_symbol:
            kZ = self.simu.get_flux_norm(machine)
        else:
            msg = "CCOutput.get_interp(): No suitable norm found."
            self.get_logger().warning(msg)
            kZ = 1

        z = [zz * kZ for zz in z]

    # interpolator = interp.interp2d(x, y, z, kind="cubic")
    interpolator = interp.CloughTocher2DInterpolator(np.array([x, y]).T, z)

    return interpolator, x, y
