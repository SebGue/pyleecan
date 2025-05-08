from ....Methods.Slot.Slot import SlotCheckError
from ....Methods.Slot.HoleM64 import *


def check(self):
    """Check that the HoleM64 object is correct

    Parameters
    ----------
    self : HoleM64
        A HoleM64 object

    Returns
    -------
    None

    Raises
    _______
    S64_WCheckError
        You must have W1+H0 <= W2
    """
    # Check that everything is set
    if self.W1 is None:
        raise S64_NoneError("You must set W1 !")
    elif self.W2 is None:
        raise S64_NoneError("You must set W2 !")
    elif self.W3 is None:
        raise S64_NoneError("You must set W3 !")
    elif self.H0 is None:
        raise S64_NoneError("You must set H0 !")
    elif self.H2 is None:
        raise S64_NoneError("You must set H2 !")
    elif self.R0 is None:
        raise S64_NoneError("You must set R0 !")

    if (2 * self.R0) >= self.H0:
        raise S64_WCheckError("You must have 2*R0 < H0")
