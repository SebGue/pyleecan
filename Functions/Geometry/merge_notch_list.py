from numpy import argsort


def merge_notch_list(notch_list_1, notch_list_2):
    """Merge and sort two notches list

    Parameters
    ----------
    notch_list_1 : list
        First notch list to merge
    notch_list_2 : list
        Second notch list to merge

    Returns
    -------
    notch_list : list 
        list of dictionary with key: "begin_angle", "end_angle", "obj"
    """
    N1 = len(notch_list_1)
    N2 = len(notch_list_2)

    # sort individual lists first
    list_1_angle = [notch["begin_angle"] for notch in notch_list_1]
    list_2_angle = [notch["begin_angle"] for notch in notch_list_2]

    sorted_id_1 = argsort(list_1_angle)
    sorted_id_2 = argsort(list_2_angle)

    list_1 = [notch_list_1[id] for id in sorted_id_1]
    list_2 = [notch_list_2[id] for id in sorted_id_2]

    merged = []
    ii, jj = 0, 0  # Index to go thought the lists

    while ii < N1 and jj < N2:
        if (
            list_1[ii]["begin_angle"] < list_2[jj]["begin_angle"]
            and list_1[ii]["end_angle"] <= list_2[jj]["begin_angle"]
        ):  # Add a notch from notch_list_1
            merged.append(list_1[ii])
            ii += 1
        elif (
            list_2[jj]["begin_angle"] < list_1[ii]["begin_angle"]
            and list_2[jj]["end_angle"] <= list_1[ii]["begin_angle"]
        ):  # Add a notch from notch_list_2
            merged.append(list_2[jj])
            jj += 1
        else:
            raise NotchError(
                "Notches and/or Slots are coliding:\n"
                + str(list_1[ii])
                + "\n"
                + str(list_2[ii])
            )

    # One of the list is not "finished"
    merged = merged + list_1[ii:] + list_2[jj:]

    return merged


class NotchError(Exception):
    """Raised when notch are coliding
    """

    pass
