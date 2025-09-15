# -*- coding: utf-8 -*-
from ...Functions.labels import BOUNDARY_PROP_LAB


def get_boundary_condition(line, boundary_prop):
    """Returns the boundary name on a line that is used in FEA coupling

    Parameters
    ----------
    line :
        a line with a label

    Returns
    -------
    label : string
        boundary name
    """
    bc_name = ""
    if line.prop_dict and BOUNDARY_PROP_LAB in line.prop_dict:
        label = line.prop_dict[BOUNDARY_PROP_LAB]
        parts = label.rsplit("-", 1)  # split from right side into 2 parts
        if len(parts) == 2 and parts[1].isdigit():
            label_without_id = parts[0]
            line_id = parts[1]
        else:
            label_without_id = label
            line_id = None
        
        if label_without_id in boundary_prop:
            bc_name = boundary_prop[label_without_id]
        
        if line_id:
            bc_name += f"_{line_id}"
    
    return bc_name
