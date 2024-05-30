# -*- coding: utf-8 -*-
from numpy import exp, pi

from ...Classes.Arc1 import Arc1
from ...Classes.Arc2 import Arc2
from ...Classes.Circle import Circle
from ...Classes.Segment import Segment
from ...Classes.SurfLine import SurfLine
from ...Functions.labels import (
    AIRBOX_R_LAB,
    NO_LAM_LAB,
    AIRBOX_LAB,
    BOUNDARY_PROP_LAB,
    AIRBOX_SR_LAB,
    AIRBOX_SL_LAB,
    AIRBOX_R_LAB,
)


def get_air_box(sym, machine):
    """Returns an outer surface surrounding the external lamination

    Parameters
    ----------
    sym: int
        Symmetry factor (1= full machine, 2= half of the machine...)
    machine: float
        to extract laminations

    Returns
    -------
    surf_list: list
        Outer surface (Air Box)
    """
    lam_list = machine.get_lam_list()
    lam_int = lam_list[0]
    lam_ext = lam_list[1]

    R_ext = lam_ext.get_Ryoke()
    R_ab = 1.5 * R_ext  # Default at 1.5 times external lam radius

    surf_list = list()
    if sym == 1:  # Complete machine
        _get_ring_segment(R_ab, R_ext, 0, pi, surf_list, idx=1)
        _get_ring_segment(R_ab, R_ext, pi, 2 * pi, surf_list, idx=2)

    else:  # Symmetry
        # Internal AirGap
        _get_ring_segment(
            R_ab,
            R_ext,
            0,
            2 * pi / sym,
            surf_list,
            prop_R={BOUNDARY_PROP_LAB: AIRBOX_SR_LAB},
            prop_L={BOUNDARY_PROP_LAB: AIRBOX_SL_LAB},
        )

    return surf_list


def _get_ring_segment(
    Rext, Rint, angle_start, angle_end, surf_list, prop_R=None, prop_L=None, idx=None
):
    """Get a ring surface."""
    Z1 = Rint * exp(1j * angle_start)
    Z0 = Rint * exp(1j * angle_end)
    Z2 = Rext * exp(1j * angle_start)
    Z3 = Rext * exp(1j * angle_end)
    airbox_lines = list()
    airbox_lines.append(Segment(begin=Z1, end=Z2, prop_dict=prop_R))
    airbox_lines.append(
        Arc1(
            begin=Z2,
            end=Z3,
            radius=Rext,
            prop_dict={BOUNDARY_PROP_LAB: AIRBOX_R_LAB},
            is_trigo_direction=True,
        )
    )
    airbox_lines.append(Segment(begin=Z3, end=Z0, prop_dict=prop_L))
    # airbox_lines.append(Arc2(begin=Z0, center=0.0, angle=-2 * pi / sym))
    airbox_lines.append(
        Arc1(
            begin=Z0,
            end=Z1,
            radius=-Rint,
            is_trigo_direction=False,
        )
    )
    lab_idx = "" if idx is None else "Pt" + str(idx)
    surf_list.append(
        SurfLine(
            line_list=airbox_lines,
            # point_ref=0.0 * Z2 * exp(1j * pi / sym),
            label=NO_LAM_LAB + "_" + AIRBOX_LAB + lab_idx,
        )
    )
