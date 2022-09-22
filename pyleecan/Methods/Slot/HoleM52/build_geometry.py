# -*- coding: utf-8 -*-

from numpy import exp

from ....Classes.Segment import Segment
from ....Classes.SurfLine import SurfLine
from ....Functions.labels import HOLEV_LAB, HOLEM_LAB
from ....Functions.labels import HOLE_LAB, HOLEM_LAB, MAG_LAB, ROTOR_LAB
from ....Functions.labels import RIGHT_LAB, LEFT_LAB, TOP_LAB, BOT_LAB
from ....Functions.labels import BOUNDARY_PROP_LAB as BND


def build_geometry(self, alpha=0, delta=0, is_simplified=False):
    """Compute the curve (Line) needed to plot the Slot.
    The ending point of a curve is the starting point of the next curve in
    the list

    Parameters
    ----------
    self : HoleM52
        A HoleM52 object
    alpha : float
        Angle to rotate the slot (Default value = 0) [rad]
    delta : complex
        Complex to translate the slot (Default value = 0)
    is_simplified : bool
       True to avoid line superposition

    Returns
    -------
    surf_list: list
        List of SurfLine needed to draw the HoleM51

    """
    # Get correct label for surfaces
    lam_label = self.parent.get_label()
    R_id, surf_type = self.get_R_id()
    vent_label = lam_label + "_" + surf_type + "_R" + str(R_id) + "-"
    mag_label = lam_label + "_" + HOLEM_LAB + "_R" + str(R_id) + "-"

    MAG = ROTOR_LAB + "_" + MAG_LAB
    HOLE = ROTOR_LAB + "_" + HOLE_LAB

    point_dict = self._comp_point_coordinate()
    Z1 = point_dict["Z1"]
    Z2 = point_dict["Z2"]
    Z3 = point_dict["Z3"]
    Z4 = point_dict["Z4"]
    Z6 = point_dict["Z6"]
    Z7 = point_dict["Z7"]
    Z8 = point_dict["Z8"]
    Z9 = point_dict["Z9"]
    Z10 = point_dict["Z10"]
    Z11 = point_dict["Z11"]

    # Creation of the air curve
    curve_list = list()
    curve_list.append(Segment(Z1, Z2))
    curve_list.append(Segment(Z2, Z3))
    curve_list.append(Segment(Z3, Z11, prop_dict={BND: "_".join([HOLE, LEFT_LAB])}))
    curve_list.append(Segment(Z11, Z1, prop_dict={BND: "_".join([HOLE, TOP_LAB, "0"])}))
    point_ref = (Z1 + Z2 + Z3 + Z11) / 4
    S1 = SurfLine(line_list=curve_list, point_ref=point_ref)

    # Creation of the magnet curve
    curve_list = list()
    if is_simplified:
        curve_list.append(Segment(Z3, Z11, prop_dict={BND: "_".join([MAG, LEFT_LAB])}))
        curve_list.append(Segment(Z7, Z10, prop_dict={BND: "_".join([MAG, RIGHT_LAB])}))
    else:
        curve_list.append(Segment(Z4, Z11, prop_dict={BND: "_".join([MAG, LEFT_LAB])}))
        curve_list.append(Segment(Z11, Z10, prop_dict={BND: "_".join([MAG, TOP_LAB])}))
        curve_list.append(Segment(Z10, Z6, prop_dict={BND: "_".join([MAG, RIGHT_LAB])}))
        curve_list.append(Segment(Z6, Z4, prop_dict={BND: "_".join([MAG, BOT_LAB])}))
    point_ref = (Z11 + Z4 + Z6 + Z10) / 4
    S2 = SurfLine(line_list=curve_list, label=mag_label + "T0-S0", point_ref=point_ref)

    # Creation of the second air curve
    curve_list = list()
    curve_list.append(Segment(Z7, Z8))
    curve_list.append(Segment(Z8, Z9))
    curve_list.append(Segment(Z9, Z10, prop_dict={BND: "_".join([HOLE, TOP_LAB, "1"])}))
    curve_list.append(Segment(Z10, Z7, prop_dict={BND: "_".join([HOLE, RIGHT_LAB])}))
    point_ref = (Z7 + Z8 + Z9 + Z10) / 4
    S3 = SurfLine(line_list=curve_list, point_ref=point_ref)

    # Area with no magnet (S1 + S2 + S3)
    curve_list = list()
    curve_list.append(Segment(Z1, Z2))
    curve_list.append(Segment(Z2, Z3))
    if self.H2 > 0:
        curve_list.append(Segment(Z3, Z4))
    curve_list.append(Segment(Z4, Z6))
    if self.H2 > 0:
        curve_list.append(Segment(Z6, Z7))
    curve_list.append(Segment(Z7, Z8))
    curve_list.append(Segment(Z8, Z9))
    curve_list.append(Segment(Z9, Z1))
    point_ref = (Z11 + Z4 + Z6 + Z10) / 4
    S4 = SurfLine(line_list=curve_list, point_ref=point_ref)

    if self.magnet_0:
        S1.label = vent_label + "T0-S0"  # Hole
        S3.label = vent_label + "T1-S0"  # Hole
        surf_list = [S1, S2, S3]
    else:
        S4.label = vent_label + "T0-S0"  # Hole
        surf_list = [S4]

    # Apply the transformations
    for surf in surf_list:
        surf.rotate(alpha)
        surf.translate(delta)

    return surf_list
