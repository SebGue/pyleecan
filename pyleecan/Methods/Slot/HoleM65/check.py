from ....Methods.Slot.Slot import SlotCheckError
from ....Methods.Slot.HoleM65 import *


def check(self):
    """Check that the HoleM65 object is correct

    Parameters
    ----------
    self : HoleM65
        A HoleM65 object

    Returns
    -------
    None

    Raises
    _______
    S65_WCheckError
        You must have W4 > R
    """
    # Check that everything is set
    if self.W1 is None:
        raise S65_NoneError("You must set W1 !")
    elif self.W2 is None:
        raise S65_NoneError("You must set W2 !")
    elif self.W3 is None:
        raise S65_NoneError("You must set W3 !")
    elif self.W4 is None:
        raise S65_NoneError("You must set W4 !")
    elif self.H0 is None:
        raise S65_NoneError("You must set H0 !")
    elif self.H1 is None:
        raise S65_NoneError("You must set H1 !")
    elif self.R is None:
        raise S65_NoneError("You must set R !")

    if self.W4 <= self.R:
        raise S65_WCheckError("You must have W4 > R")

    if self.H1 <= (2 * self.R):
        raise S65_WCheckError("You must have H1 > 2 * R")
