from numpy import pi, exp, angle

from ....Functions.Geometry.cut_lines_between_angle import cut_lines_between_angles
from ....Classes.Segment import Segment


def merge_slot_translate(self, radius_desc_list, prop_dict, sym):
    """Merge the Bore shape with notches/slot on the bore/yoke
    Translate method: Translate lines of notch/slot to match the radius of the shape
    technically Keep the slot/notch height by reduicing the yoke height

    Parameters
    ----------
    radius_desc_list : list
        List of dict to describe the bore/yoke radius
    prop_dict : dict
        Property dictionary to apply on the radius lines (not on slot/notch)
    sym : int
        Symmetry factor (1= full machine, 2= half of the machine...)

    Returns
    -------
    line_list : list
        List of lines needed to draw the radius
    """
    # Get the actucal bore lines (0 to 2*pi) without nocth/slot
    radius_lines = self.get_bore_line()

    # Update begin and end angle if Radius (next step already cut lines)
    if sym != 1 and radius_desc_list[0]["label"] == "Radius":
        radius_desc_list[0]["begin_angle"] = 0
    if sym != 1 and radius_desc_list[-1]["label"] == "Radius":
        radius_desc_list[-1]["end_angle"] = 2 * pi / sym

    # Apply merge strategy on slot/notch
    line_list = list()

    # first translate the slot/notch lines to match the bore
    for desc_dict in radius_desc_list:
        if desc_dict["label"] != "Radius":  # notch/slot
            lines = desc_dict["lines"]
            # create 2 lines parallel to the slot radius line on the opening points
            begin = lines[0].get_begin()
            end = lines[-1].get_end()
            mid = (begin + end) / 2
            ang = angle(mid)
            line1 = Segment(begin=begin - mid, end=begin + mid)
            line2 = Segment(begin=end - mid, end=end + mid)

            # find intersection points of lines with bore
            inter1, inter2 = list(), list()
            for line in radius_lines:
                Zi1 = line.intersect_line(line1.get_begin(), line1.get_end())
                Zi2 = line.intersect_line(line2.get_begin(), line2.get_end())
                inter1.extend([z for z in Zi1 if (z * exp(-1j * ang)).real > 0])
                inter2.extend([z for z in Zi2 if (z * exp(-1j * ang)).real > 0])

            if len(inter1) != 1 or len(inter2) != 1:
                raise ValueError("Error placing Notch on Bore")  # TODO: proper err.

            # do the translation of the slot lines
            dZ1 = inter1[0] - begin
            dZ2 = inter2[0] - end
            dZ = dZ1 if abs(dZ1) > abs(dZ2) else dZ2

            for line in lines:
                line.translate(dZ)

            # add line to complete the slot
            if abs(dZ1) > abs(dZ2):
                lines.append(Segment(begin=lines[-1].get_end(), end=inter2[0]))
            elif abs(dZ1) < abs(dZ2):
                lines.insert(0, Segment(begin=inter1[0], end=lines[0].get_begin()))

            # update angles
            desc_dict["begin_angle"] = angle(lines[0].get_begin())
            desc_dict["end_angle"] = angle(lines[-1].get_end())

    # Replace Arc radius from desc by lines from actual bore shape
    for ii, desc_dict in enumerate(radius_desc_list):
        if desc_dict["label"] == "Radius":
            if ii != 0 and radius_desc_list[ii - 1]["label"] != "Radius":
                begin_ang = angle(radius_desc_list[ii - 1]["lines"][-1].get_end())
            else:
                begin_ang = desc_dict["begin_angle"]

            if ii == len(radius_desc_list) - 1:
                end_ang = desc_dict["end_angle"]
            elif radius_desc_list[ii + 1]["label"] != "Radius":
                end_ang = angle(radius_desc_list[ii + 1]["lines"][0].get_begin())

            end_ang = desc_dict["end_angle"]
            desc_dict["lines"] = cut_lines_between_angles(
                radius_lines, begin_ang, end_ang
            )
            # Add prop_dict on all the Radius Lines
            if prop_dict is not None:
                for line in desc_dict["lines"]:
                    if line.prop_dict is None:
                        line.prop_dict = dict()
                    line.prop_dict.update(prop_dict)

            # update angles
            desc_dict["begin_angle"] = angle(lines[0].get_begin())
            desc_dict["end_angle"] = angle(lines[-1].get_end())

    # If slot/notch are coliding with sym lines => Cut
    if sym != 1:
        # Cut first desc (if needed)
        if radius_desc_list[0]["begin_angle"] < 0:
            lines = list()
            for line in radius_desc_list[0]["lines"]:
                top_split_list, _ = line.split_line(0, 1)
                lines.extend(top_split_list)
            radius_desc_list[0]["begin_angle"] = 0
            radius_desc_list[0]["lines"] = lines
        # Cut last desc (if needed)
        if radius_desc_list[-1]["end_angle"] > 2 * pi / sym:
            lines = list()
            for line in radius_desc_list[-1]["lines"]:
                _, bot_split_list = line.split_line(0, exp(1j * 2 * pi / sym))
                lines.extend(bot_split_list)
            radius_desc_list[-1]["end_angle"] = 2 * pi / sym
            radius_desc_list[-1]["lines"] = lines

    # create the line_list
    for desc_dict in radius_desc_list:
        line_list.extend(desc_dict["lines"])

    return line_list

    """ org. code

    notch_list = list()
    op = self.notch_shape.comp_angle_opening()

    notch = self.notch_shape.build_geometry()
    mid = (notch[0].get_begin() + notch[-1].get_end()) / 2
    line1 = Segment(begin=notch[0].get_begin() - mid, end=notch[0].get_begin())
    line2 = Segment(begin=notch[-1].get_end() - mid, end=notch[-1].get_end())

    bore_lines = self.parent.get_bore_line()
    ang = 2 * pi / self.notch_shape.Zs

    for ii in range(self.notch_shape.Zs // sym):
        ang_rot = ang * ii + self.alpha
        # find intersection points of lines with bore
        inter1, inter2 = list(), list()
        for line in bore_lines:
            line = line.copy()
            line.rotate(-ang_rot)
            Zi1 = line.intersect_line(line1.get_begin(), line1.get_end())
            Zi2 = line.intersect_line(line2.get_begin(), line2.get_end())
            inter1.extend([z for z in Zi1 if z.real > 0])
            inter2.extend([z for z in Zi2 if z.real > 0])

        if len(inter1) != 1 or len(inter2) != 1:
            raise ValueError("Error placing Notch on Bore")  # TODO: proper err.

        # translate notch
        dZ1 = inter1[0] - line1.end
        dZ2 = inter2[0] - line2.end
        dZ = dZ1 if abs(dZ1) > abs(dZ2) else dZ2

        _notch = list()
        for line in notch:
            line = line.copy()
            line.translate(dZ)
            line.rotate(ang_rot)
            _notch.append(line)

        # add line to complete notch
        if abs(dZ1) > abs(dZ2):
            _notch.append(
                Segment(begin=_notch[-1].get_end(), end=inter2[0] * exp(1j * ang_rot))
            )
        elif abs(dZ1) < abs(dZ2):
            _notch.insert(
                0,
                Segment(begin=inter1[0] * exp(1j * ang_rot), end=_notch[0].get_begin()),
            )

        # add notch to list
        notch_dict = dict()
        notch_dict["begin_angle"] = angle(_notch[0].get_begin())
        notch_dict["end_angle"] = angle(_notch[-1].get_end())
        notch_dict["obj"] = _notch
        notch_list.append(notch_dict)

    """
