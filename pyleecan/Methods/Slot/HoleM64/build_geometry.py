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

    Z1 = point_dict["Z1"]
    Z21 = point_dict["Z21"]
    Z22 = point_dict["Z22"]
    Z2c = point_dict["Z2c"]
    Z31 = point_dict["Z31"]
    Z32 = point_dict["Z32"]
    Z3c = point_dict["Z3c"]
    Z4 = point_dict["Z4"]

    ZM1 = point_dict["ZM1"]
    ZM2 = point_dict["ZM2"]
    ZM3 = point_dict["ZM3"]
    ZM4 = point_dict["ZM4"]

    # Create all the surfaces for all the cases

    # upper surface
    curve_list = list()
    curve_list.append(Segment(Z1, Z21))
    if self.R0 != 0:
        curve_list.append(Arc2(begin=Z21, center=Z2c, angle=pi / 2))

    if self.W1 > self.R0:
        curve_list.append(Segment(Z22, ZM2))
    curve_list.append(Segment(ZM2, ZM1))
    curve_list.append(Segment(ZM1, Z1))

    point_ref = (Z1 + ZM2) / 2
    S0 = SurfLine(line_list=curve_list, point_ref=point_ref)
    S0.label = vent_label + "T0-S0"

    S0s = _mirror_rotate(S0, self.Zh, vent_label + "T2-S0")

    # lower surface
    curve_list = list()
    if self.W3 > self.R0:
        curve_list.append(Segment(ZM3, Z31))

    if self.R0 != 0:
        curve_list.append(Arc2(begin=Z31, center=Z3c, angle=pi / 2))
    curve_list.append(Segment(Z32, Z4))
    curve_list.append(Segment(Z4, ZM4))
    curve_list.append(Segment(ZM4, ZM3))
    point_ref = (Z4 + ZM3) / 2
    S1 = SurfLine(line_list=curve_list, point_ref=point_ref)
    S1.label = vent_label + "T1-S0"

    S1s = _mirror_rotate(S1, self.Zh, vent_label + "T3-S0")

    # magnet_0 surface
    if self.magnet_0:
        curve_list = list()
        curve_list.append(Segment(ZM1, ZM2))
        curve_list.append(Segment(ZM2, ZM3))
        curve_list.append(Segment(ZM3, ZM4))
        curve_list.append(Segment(ZM4, ZM1))
        point_ref = (ZM1 + ZM3) / 2
        SM = SurfLine(line_list=curve_list, point_ref=point_ref)
        SM.label = mag_label + "T0-S0"

        SMs = _mirror_rotate(SM, self.Zh, mag_label + "T1-S0")

    # Create the surface list by selecting the correct ones
    surf_list = list()
    surf_list += [S0]
    if self.magnet_0:
        surf_list += [SM]
    surf_list += [S1]

    surf_list += [S0s]
    if self.magnet_0:
        surf_list += [SMs]
    surf_list += [S1s]

    # Apply the transformations
    for surf in surf_list:
        surf.rotate(alpha - pi / self.Zh)
        surf.translate(delta)

    return surf_list


def _mirror_rotate(surf, Zh, label):
    s = surf.copy()
    s.label = label
    rotate = exp(1j * 2 * pi / Zh)
    for line in s.line_list:
        line.begin = conjugate(line.begin) * rotate
        if isinstance(line, Arc2):
            line.center = conjugate(line.center) * rotate
            line.angle = -line.angle
        else:
            line.end = conjugate(line.end) * rotate

    s.point_ref = conjugate(s.point_ref) * rotate
    return s
