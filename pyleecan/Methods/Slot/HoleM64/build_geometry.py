from numpy import pi, conjugate, exp

from pyleecan.Classes.Arc2 import Arc2
from pyleecan.Classes.Segment import Segment
from pyleecan.Classes.SurfLine import SurfLine
from pyleecan.Functions.labels import HOLEM_LAB


def build_geometry(self, alpha=0, delta=0, is_simplified=False):
    """Compute the curve (Segment) needed to plot the Slot.
    The ending point of a curve is the starting point of the next curve in
    the list

    Parameters
    ----------
    self : HoleM64
        A HoleM64 object
    alpha : float
        Angle to rotate the slot (Default value = 0) [rad]
    delta : complex
        Complex to translate the slot (Default value = 0)
    is_simplified : bool
       True to avoid line superposition

    Returns
    -------
    surf_list: list
        List of SurfLine needed to draw the HoleM58

    """

    # Get correct label for surfaces
    lam_label = self.parent.get_label()
    R_id, surf_type = self.get_R_id()
    vent_label = lam_label + "_" + surf_type + "_R" + str(R_id) + "-"
    mag_label = lam_label + "_" + HOLEM_LAB + "_R" + str(R_id) + "-"

    # Get all the points
    point_dict = self._comp_point_coordinate()

    Z11 = point_dict["Z11"]
    Z12 = point_dict["Z12"]
    Z1c = point_dict["Z1c"]

    Z21 = point_dict["Z21"]
    Z22 = point_dict["Z22"]
    Z2c = point_dict["Z2c"]

    ZM1 = point_dict["ZM1"]
    ZM2 = point_dict["ZM2"]

    # Create all the surfaces for all the cases

    # upper surface
    upper_lines = list()
    if self.W1 > self.R0:
        upper_lines.append(Segment(conjugate(ZM1), conjugate(Z12)))
    if self.R0 != 0:
        upper_lines.append(
            Arc2(begin=conjugate(Z12), center=conjugate(Z1c), angle=pi / 2)
        )

    upper_lines.append(Segment(conjugate(Z11), Z11))
    if self.R0 != 0:
        upper_lines.append(Arc2(begin=Z11, center=Z1c, angle=pi / 2))

    if self.W1 > self.R0:
        upper_lines.append(Segment(Z12, ZM1))
    upper_lines.append(Segment(ZM1, conjugate(ZM1)))

    point_ref = (Z11 + conjugate(ZM1)) / 2
    S0 = SurfLine(line_list=upper_lines, point_ref=point_ref)
    S0.label = vent_label + "T0-S0"

    # lower surface
    lower_lines = list()
    if self.W3 > self.R0:
        lower_lines.append(Segment(ZM2, Z21))

    if self.R0 != 0:
        lower_lines.append(Arc2(begin=Z21, center=Z2c, angle=pi / 2))
    lower_lines.append(Segment(Z22, conjugate(Z22)))
    if self.R0 != 0:
        lower_lines.append(
            Arc2(begin=conjugate(Z22), center=conjugate(Z2c), angle=pi / 2)
        )
    if self.W3 > self.R0:
        lower_lines.append(Segment(conjugate(Z21), conjugate(ZM2)))
    lower_lines.append(Segment(conjugate(ZM2), ZM2))

    point_ref = (Z21 + conjugate(ZM2)) / 2
    S1 = SurfLine(line_list=lower_lines, point_ref=point_ref)
    S1.label = vent_label + "T1-S0"

    # magnet_0 surface
    magnet_lines = list()
    magnet_lines.append(Segment(conjugate(ZM1), ZM1))
    magnet_lines.append(Segment(ZM1, ZM2))
    magnet_lines.append(Segment(ZM2, conjugate(ZM2)))
    magnet_lines.append(Segment(conjugate(ZM2), conjugate(ZM1)))
    point_ref = (ZM1 + conjugate(ZM2)) / 2
    SM = SurfLine(line_list=magnet_lines, point_ref=point_ref)
    SM.label = mag_label + "T0-S0"

    # hole without magnet_0
    hole_lines = upper_lines[:-2]
    hole_lines += [magnet_lines[1]]
    hole_lines += lower_lines[:-2]
    hole_lines += [magnet_lines[3]]
    point_ref = (ZM1 + conjugate(ZM2)) / 2
    S10 = SurfLine(line_list=hole_lines, point_ref=point_ref)
    S10.label = vent_label + "T0-S0"

    # Create the surface list by selecting the correct ones
    if self.magnet_0:
        surf_list = [S0, SM, S1]
    else:
        surf_list = [S10]

    # Apply the transformations
    for surf in surf_list:
        surf.rotate(alpha - 1 * pi / self.Zh)
        surf.translate(delta)

    return surf_list
