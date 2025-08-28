from numpy import pi
from ....Classes.Segment import Segment


def comp_magnetization_dict(self, is_north=True):
    """Compute the dictionary of the magnetization direction of the magnets (key=magnet_X, value=angle[rad])
    Mangetization angle with Hole centered on Ox axis

    Parameters
    ----------
    self : HoleM64
        a HoleM64 object
    is_north: True
        True: comp north magnetization, else add pi [rad]

    Returns
    -------
    mag_dict: dict
        magnetization dictionary (key=magnet_X, value=angle[rad])
    """

    # Comp magnet
    point_dict = self._comp_point_coordinate()

    mag_dict = dict()
    S0 = Segment(point_dict["ZM1"], point_dict["ZM2"])
    mag_dict["magnet_0"] = S0.comp_normal() - pi / self.Zh + pi
    # mag_dict["magnet_1"] = mag_dict["magnet_0"] + 2 * pi / self.Zh + pi

    if not is_north:
        mag_dict["magnet_0"] += 2 * pi / self.Zh + pi
        # mag_dict["magnet_1"] += pi

    if self.magnetization_dict_offset is not None:
        for key, value in self.magnetization_dict_offset:
            mag_dict[key] += value

    return mag_dict
