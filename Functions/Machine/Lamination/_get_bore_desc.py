from pyleecan.Functions.Geometry.merge_notch_list import merge_notch_list
from pyleecan.Classes.Arc1 import Arc1
from numpy import exp, pi


def _get_bore_desc(lam, sym=1):
    """This function returns an ordered description of the elements 
    that defines the bore radius of the lamination

    Parameters
    ----------
    lam : Lamination
        A Lamination object
    
    sym : int
        Symmetry factor (1= full machine, 2= half of the machine...)

    Returns
    -------
    bore_desc : list 
        list of dictionary with key: "begin_angle", "end_angle", "obj"
    """
    Rbo = lam.get_Rbo()

    # get the slots
    if hasattr(lam, "slot_list") and hasattr(lam, "alpha"):
        slots = lam.slot_list
        alphas = lam.alpha
    elif hasattr(lam, "slot"):
        slots = [lam.slot]
        alphas = [pi / lam.slot.Zs]
    else:
        slots = []
        alphas = []

    slot_list = list()
    # First add all the slots
    for slot, alpha in zip(slots, alphas):
        op = slot.comp_angle_opening()
        for ii in range(slot.Zs // sym):
            bore_dict = dict()
            beta = 2 * pi * ii / slot.Zs
            bore_dict["begin_angle"] = alpha - op / 2 + beta
            bore_dict["end_angle"] = alpha + op / 2 + beta
            bore_dict["obj"] = slot
            slot_list.append(bore_dict)

    # Get the notches
    notch_list = lam.get_notch_list(sym=sym)

    # Merge Slot and Notches
    merged_list = merge_notch_list(slot_list, notch_list)

    # Add all the bore lines
    bore_desc = list()
    # if lamination has slots and/or notches
    if merged_list:
        for ii, desc in enumerate(merged_list):
            bore_desc.append(desc)
            if ii != len(merged_list) - 1:
                bore_dict = dict()
                bore_dict["begin_angle"] = merged_list[ii]["end_angle"]
                bore_dict["end_angle"] = merged_list[ii + 1]["begin_angle"]
                bore_dict["obj"] = Arc1(
                    begin=Rbo * exp(1j * bore_dict["begin_angle"]),
                    end=Rbo * exp(1j * bore_dict["end_angle"]),
                    radius=Rbo,
                    is_trigo_direction=True,
                )
                bore_desc.append(bore_dict)
    
        # Add last bore line
        if sym == 1:
            bore_dict = dict()
            bore_dict["begin_angle"] = merged_list[-1]["end_angle"]
            bore_dict["end_angle"] = merged_list[0]["begin_angle"]
            bore_dict["obj"] = Arc1(
                begin=Rbo * exp(1j * bore_dict["begin_angle"]),
                end=Rbo * exp(1j * bore_dict["end_angle"]),
                radius=Rbo,
                is_trigo_direction=True,
            )
            if merged_list[0]["begin_angle"] < 0:
                # First element is an slot or notch
                bore_desc.append(bore_dict)
            else:
                # First element is a bore line
                bore_desc.insert(0, bore_dict)
        else:  # With symmetry
            # Add last bore line
            bore_dict = dict()
            bore_dict["begin_angle"] = merged_list[-1]["end_angle"]
            bore_dict["end_angle"] = 2 * pi / sym
            bore_dict["obj"] = Arc1(
                begin=Rbo * exp(1j * bore_dict["begin_angle"]),
                end=Rbo * exp(1j * bore_dict["end_angle"]),
                radius=Rbo,
                is_trigo_direction=True,
            )
            bore_desc.append(bore_dict)

            # Add first bore line
            bore_dict = dict()
            bore_dict["begin_angle"] = 0
            bore_dict["end_angle"] = merged_list[0]["begin_angle"]
            bore_dict["obj"] = Arc1(
                begin=Rbo * exp(1j * bore_dict["begin_angle"]),
                end=Rbo * exp(1j * bore_dict["end_angle"]),
                radius=Rbo,
                is_trigo_direction=True,
            )
            bore_desc.insert(0, bore_dict)
    
    # if lamination has no notches or slots
    else:
        if sym == 1:
            # first half circle
            bore_dict = dict()
            bore_dict["begin_angle"] = 0
            bore_dict["end_angle"] = pi
            bore_dict["obj"] = Arc1(
                    begin=Rbo * exp(1j * bore_dict["begin_angle"]),
                    end=Rbo * exp(1j * bore_dict["end_angle"]),
                    radius=Rbo,
                    is_trigo_direction=True,
                )
            bore_desc.append(bore_dict)
            # second half circle
            bore_dict = dict()
            bore_dict["begin_angle"] = pi
            bore_dict["end_angle"] = 2 * pi
            bore_dict["obj"] = Arc1(
                    begin=Rbo * exp(1j * bore_dict["begin_angle"]),
                    end=Rbo * exp(1j * bore_dict["end_angle"]),
                    radius=Rbo,
                    is_trigo_direction=True,
                )
            bore_desc.append(bore_dict)
        else:
            bore_dict = dict()
            bore_dict["begin_angle"] = 0
            bore_dict["end_angle"] = 2 * pi / sym
            bore_dict["obj"] = Arc1(
                    begin=Rbo * exp(1j * bore_dict["begin_angle"]),
                    end=Rbo * exp(1j * bore_dict["end_angle"]),
                    radius=Rbo,
                    is_trigo_direction=True,
                )
            bore_desc.append(bore_dict)
            

    return bore_desc
