# -*- coding: utf-8 -*-
# File generated according to Generator/ClassesRef/Simulation/Electrical.csv
# WARNING! All changes made in this file will be lost!
"""Method code available at https://github.com/Eomys/pyleecan/tree/master/pyleecan/Methods/Simulation/Electrical
"""

from os import linesep
from logging import getLogger
from ._check import check_var, raise_
from ..Functions.get_logger import get_logger
from ..Functions.save import save
from ..Functions.copy import copy
from ..Functions.load import load_init_dict
from ..Functions.Load.import_class import import_class
from ._frozen import FrozenClass

# Import all class method
# Try/catch to remove unnecessary dependencies in unused method
try:
    from ..Methods.Simulation.Electrical.run import run
except ImportError as error:
    run = error

try:
    from ..Methods.Simulation.Electrical.comp_power import comp_power
except ImportError as error:
    comp_power = error

try:
    from ..Methods.Simulation.Electrical.comp_torque import comp_torque
except ImportError as error:
    comp_torque = error


from ._check import InitUnKnowClassError
from .EEC import EEC


class Electrical(FrozenClass):
    """Electric module object for electrical equivalent circuit simulation"""

    VERSION = 1

    # Check ImportError to remove unnecessary dependencies in unused method
    # cf Methods.Simulation.Electrical.run
    if isinstance(run, ImportError):
        run = property(
            fget=lambda x: raise_(
                ImportError("Can't use Electrical method run: " + str(run))
            )
        )
    else:
        run = run
    # cf Methods.Simulation.Electrical.comp_power
    if isinstance(comp_power, ImportError):
        comp_power = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use Electrical method comp_power: " + str(comp_power)
                )
            )
        )
    else:
        comp_power = comp_power
    # cf Methods.Simulation.Electrical.comp_torque
    if isinstance(comp_torque, ImportError):
        comp_torque = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use Electrical method comp_torque: " + str(comp_torque)
                )
            )
        )
    else:
        comp_torque = comp_torque
    # save and copy methods are available in all object
    save = save
    copy = copy
    # get_logger method is available in all object
    get_logger = get_logger

    def __init__(
        self,
        eec=None,
        is_comp_torque=False,
        Imax=None,
        Umax=None,
        init_dict=None,
        init_str=None,
    ):
        """Constructor of the class. Can be use in three ways :
        - __init__ (arg1 = 1, arg3 = 5) every parameters have name and default values
            for pyleecan type, -1 will call the default constructor
        - __init__ (init_dict = d) d must be a dictionnary with property names as keys
        - __init__ (init_str = s) s must be a string
        s is the file path to load

        ndarray or list can be given for Vector and Matrix
        object or dict can be given for pyleecan Object"""

        if init_str is not None:  # Load from a file
            init_dict = load_init_dict(init_str)[1]
        if init_dict is not None:  # Initialisation by dict
            assert type(init_dict) is dict
            # Overwrite default value with init_dict content
            if "eec" in list(init_dict.keys()):
                eec = init_dict["eec"]
            if "is_comp_torque" in list(init_dict.keys()):
                is_comp_torque = init_dict["is_comp_torque"]
            if "Imax" in list(init_dict.keys()):
                Imax = init_dict["Imax"]
            if "Umax" in list(init_dict.keys()):
                Umax = init_dict["Umax"]
        # Set the properties (value check and convertion are done in setter)
        self.parent = None
        self.eec = eec
        self.is_comp_torque = is_comp_torque
        self.Imax = Imax
        self.Umax = Umax

        # The class is frozen, for now it's impossible to add new properties
        self._freeze()

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        Electrical_str = ""
        if self.parent is None:
            Electrical_str += "parent = None " + linesep
        else:
            Electrical_str += "parent = " + str(type(self.parent)) + " object" + linesep
        if self.eec is not None:
            tmp = self.eec.__str__().replace(linesep, linesep + "\t").rstrip("\t")
            Electrical_str += "eec = " + tmp
        else:
            Electrical_str += "eec = None" + linesep + linesep
        Electrical_str += "is_comp_torque = " + str(self.is_comp_torque) + linesep
        Electrical_str += "Imax = " + str(self.Imax) + linesep
        Electrical_str += "Umax = " + str(self.Umax) + linesep
        return Electrical_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False
        if other.eec != self.eec:
            return False
        if other.is_comp_torque != self.is_comp_torque:
            return False
        if other.Imax != self.Imax:
            return False
        if other.Umax != self.Umax:
            return False
        return True

    def as_dict(self):
        """Convert this object in a json seriable dict (can be use in __init__)"""

        Electrical_dict = dict()
        if self.eec is None:
            Electrical_dict["eec"] = None
        else:
            Electrical_dict["eec"] = self.eec.as_dict()
        Electrical_dict["is_comp_torque"] = self.is_comp_torque
        Electrical_dict["Imax"] = self.Imax
        Electrical_dict["Umax"] = self.Umax
        # The class name is added to the dict for deserialisation purpose
        Electrical_dict["__class__"] = "Electrical"
        return Electrical_dict

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""

        if self.eec is not None:
            self.eec._set_None()
        self.is_comp_torque = None
        self.Imax = None
        self.Umax = None

    def _get_eec(self):
        """getter of eec"""
        return self._eec

    def _set_eec(self, value):
        """setter of eec"""
        if isinstance(value, str):  # Load from file
            value = load_init_dict(value)[1]
        if isinstance(value, dict) and "__class__" in value:
            class_obj = import_class("pyleecan.Classes", value.get("__class__"), "eec")
            value = class_obj(init_dict=value)
        elif type(value) is int and value == -1:  # Default constructor
            value = EEC()
        check_var("eec", value, "EEC")
        self._eec = value

        if self._eec is not None:
            self._eec.parent = self

    eec = property(
        fget=_get_eec,
        fset=_set_eec,
        doc=u"""Electrical Equivalent Circuit

        :Type: EEC
        """,
    )

    def _get_is_comp_torque(self):
        """getter of is_comp_torque"""
        return self._is_comp_torque

    def _set_is_comp_torque(self, value):
        """setter of is_comp_torque"""
        check_var("is_comp_torque", value, "bool")
        self._is_comp_torque = value

    is_comp_torque = property(
        fget=_get_is_comp_torque,
        fset=_set_is_comp_torque,
        doc=u"""Iteratively compute currents and voltages coressponding to reference torque

        :Type: bool
        """,
    )

    def _get_Imax(self):
        """getter of Imax"""
        return self._Imax

    def _set_Imax(self, value):
        """setter of Imax"""
        check_var("Imax", value, "float", Vmin=0)
        self._Imax = value

    Imax = property(
        fget=_get_Imax,
        fset=_set_Imax,
        doc=u"""Maximum phase peak current

        :Type: float
        :min: 0
        """,
    )

    def _get_Umax(self):
        """getter of Umax"""
        return self._Umax

    def _set_Umax(self, value):
        """setter of Umax"""
        check_var("Umax", value, "float", Vmin=0)
        self._Umax = value

    Umax = property(
        fget=_get_Umax,
        fset=_set_Umax,
        doc=u"""Maximum phase peak voltage

        :Type: float
        :min: 0
        """,
    )
