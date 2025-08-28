from ....Methods.Geometry.Segment import *
from ....Functions.labels import BOUNDARY_PROP_LAB


def check(self):
    """assert that the line is correct (begin != end)

    Parameters
    ----------
    self : Segment
        A Segment object

    Returns
    -------
    None

    Raises
    ------
    PointSegmentError
        The beginning point and the ending point of an Segment
        can't be the same

    """
    if self.begin == self.end or (abs(self.begin) == 0 and abs(self.end) == 0):
        info = ""
        if self.prop_dict:
            if self.parent and self.parent.label:
                info += f"Parent surface: '{self.parent.label}'\n"
            info += "Segement properties:\n"
            for key, value in self.prop_dict.items():
                info += f"{key}: {value}\n"

        raise PointSegmentError(
            "The beginning point and the ending point of a Segment can't be the same.\n"
            + info
        )
