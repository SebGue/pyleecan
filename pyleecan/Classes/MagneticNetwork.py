# -*- coding: utf-8 -*-
# File generated according to Generator/ClassesRef/Simulation/MagneticNetwork.csv
# WARNING! All changes made in this file will be lost!
"""Method code available at https://github.com/Eomys/pyleecan/tree/master/pyleecan/Methods/Simulation/MagneticNetwork
"""

from os import linesep
from sys import getsizeof
from logging import getLogger
from ._check import check_var, raise_
from ..Functions.get_logger import get_logger
from ..Functions.save import save
from ..Functions.load import load_init_dict
from ..Functions.Load.import_class import import_class
from copy import deepcopy
from ._frozen import FrozenClass

# Import all class method
# Try/catch to remove unnecessary dependencies in unused method
try:
    from ..Methods.Simulation.MagneticNetwork.assembler import assembler
except ImportError as error:
    assembler = error

try:
    from ..Methods.Simulation.MagneticNetwork.geometry_linear_motor import (
        geometry_linear_motor,
    )
except ImportError as error:
    geometry_linear_motor = error

try:
    from ..Methods.Simulation.MagneticNetwork.geometry_linear_motor_separetion import (
        geometry_linear_motor_separetion,
    )
except ImportError as error:
    geometry_linear_motor_separetion = error

try:
    from ..Methods.Simulation.MagneticNetwork.plot import plot
except ImportError as error:
    plot = error

try:
    from ..Methods.Simulation.MagneticNetwork.post_processing import post_processing
except ImportError as error:
    post_processing = error

try:
    from ..Methods.Simulation.MagneticNetwork.pre_processing import pre_processing
except ImportError as error:
    pre_processing = error

try:
    from ..Methods.Simulation.MagneticNetwork.RN_linear_motor import RN_linear_motor
except ImportError as error:
    RN_linear_motor = error

try:
    from ..Methods.Simulation.MagneticNetwork.run import run
except ImportError as error:
    run = error

try:
    from ..Methods.Simulation.MagneticNetwork.run_non_linear import run_non_linear
except ImportError as error:
    run_non_linear = error

try:
    from ..Methods.Simulation.MagneticNetwork.run_radial import run_radial
except ImportError as error:
    run_radial = error

try:
    from ..Methods.Simulation.MagneticNetwork.solver_linear_model import (
        solver_linear_model,
    )
except ImportError as error:
    solver_linear_model = error

try:
    from ..Methods.Simulation.MagneticNetwork.solver_non_linear_model import (
        solver_non_linear_model,
    )
except ImportError as error:
    solver_non_linear_model = error

try:
    from ..Methods.Simulation.MagneticNetwork.solver_plus_non_linear_model.py import py
except ImportError as error:
    py = error


from numpy import isnan
from ._check import InitUnKnowClassError


class MagneticNetwork(FrozenClass):
    """Abstract class to solve the electric machines using the reluctance network"""

    VERSION = 1

    # Check ImportError to remove unnecessary dependencies in unused method
    # cf Methods.Simulation.MagneticNetwork.assembler
    if isinstance(assembler, ImportError):
        assembler = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method assembler: " + str(assembler)
                )
            )
        )
    else:
        assembler = assembler
    # cf Methods.Simulation.MagneticNetwork.geometry_linear_motor
    if isinstance(geometry_linear_motor, ImportError):
        geometry_linear_motor = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method geometry_linear_motor: "
                    + str(geometry_linear_motor)
                )
            )
        )
    else:
        geometry_linear_motor = geometry_linear_motor
    # cf Methods.Simulation.MagneticNetwork.geometry_linear_motor_separetion
    if isinstance(geometry_linear_motor_separetion, ImportError):
        geometry_linear_motor_separetion = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method geometry_linear_motor_separetion: "
                    + str(geometry_linear_motor_separetion)
                )
            )
        )
    else:
        geometry_linear_motor_separetion = geometry_linear_motor_separetion
    # cf Methods.Simulation.MagneticNetwork.plot
    if isinstance(plot, ImportError):
        plot = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagneticNetwork method plot: " + str(plot))
            )
        )
    else:
        plot = plot
    # cf Methods.Simulation.MagneticNetwork.post_processing
    if isinstance(post_processing, ImportError):
        post_processing = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method post_processing: "
                    + str(post_processing)
                )
            )
        )
    else:
        post_processing = post_processing
    # cf Methods.Simulation.MagneticNetwork.pre_processing
    if isinstance(pre_processing, ImportError):
        pre_processing = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method pre_processing: "
                    + str(pre_processing)
                )
            )
        )
    else:
        pre_processing = pre_processing
    # cf Methods.Simulation.MagneticNetwork.RN_linear_motor
    if isinstance(RN_linear_motor, ImportError):
        RN_linear_motor = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method RN_linear_motor: "
                    + str(RN_linear_motor)
                )
            )
        )
    else:
        RN_linear_motor = RN_linear_motor
    # cf Methods.Simulation.MagneticNetwork.run
    if isinstance(run, ImportError):
        run = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagneticNetwork method run: " + str(run))
            )
        )
    else:
        run = run
    # cf Methods.Simulation.MagneticNetwork.run_non_linear
    if isinstance(run_non_linear, ImportError):
        run_non_linear = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method run_non_linear: "
                    + str(run_non_linear)
                )
            )
        )
    else:
        run_non_linear = run_non_linear
    # cf Methods.Simulation.MagneticNetwork.run_radial
    if isinstance(run_radial, ImportError):
        run_radial = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method run_radial: " + str(run_radial)
                )
            )
        )
    else:
        run_radial = run_radial
    # cf Methods.Simulation.MagneticNetwork.solver_linear_model
    if isinstance(solver_linear_model, ImportError):
        solver_linear_model = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method solver_linear_model: "
                    + str(solver_linear_model)
                )
            )
        )
    else:
        solver_linear_model = solver_linear_model
    # cf Methods.Simulation.MagneticNetwork.solver_non_linear_model
    if isinstance(solver_non_linear_model, ImportError):
        solver_non_linear_model = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagneticNetwork method solver_non_linear_model: "
                    + str(solver_non_linear_model)
                )
            )
        )
    else:
        solver_non_linear_model = solver_non_linear_model
    # cf Methods.Simulation.MagneticNetwork.solver_plus_non_linear_model.py
    if isinstance(py, ImportError):
        py = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagneticNetwork method py: " + str(py))
            )
        )
    else:
        py = py
    # generic save method is available in all object
    save = save
    # get_logger method is available in all object
    get_logger = get_logger

    def __init__(self, init_dict=None, init_str=None):
        """Constructor of the class. Can be use in two ways :
        - __init__ (arg1 = 1, arg3 = 5) every parameters have name and default values
            for Matrix, None will initialise the property with an empty Matrix
            for pyleecan type, None will call the default constructor
        - __init__ (init_dict = d) d must be a dictionary wiht every properties as keys

        ndarray or list can be given for Vector and Matrix
        object or dict can be given for pyleecan Object"""

        if init_dict is not None:  # Initialisation by dict
            assert "__class__" in init_dict
            assert init_dict["__class__"] == "MagneticNetwork"
        if init_str is not None:  # Initialisation by str
            assert type(init_str) is str
        # The class is frozen, for now it's impossible to add new properties
        self.parent = None
        self._freeze()

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        MagneticNetwork_str = ""
        if self.parent is None:
            MagneticNetwork_str += "parent = None " + linesep
        else:
            MagneticNetwork_str += (
                "parent = " + str(type(self.parent)) + " object" + linesep
            )
        return MagneticNetwork_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False
        return True

    def compare(self, other, name="self", ignore_list=None, is_add_value=False):
        """Compare two objects and return list of differences"""

        if ignore_list is None:
            ignore_list = list()
        if type(other) != type(self):
            return ["type(" + name + ")"]
        diff_list = list()
        # Filter ignore differences
        diff_list = list(filter(lambda x: x not in ignore_list, diff_list))
        return diff_list

    def __sizeof__(self):
        """Return the size in memory of the object (including all subobject)"""

        S = 0  # Full size of the object
        return S

    def as_dict(self, type_handle_ndarray=0, keep_function=False, **kwargs):
        """
        Convert this object in a json serializable dict (can be use in __init__).
        type_handle_ndarray: int
            How to handle ndarray (0: tolist, 1: copy, 2: nothing)
        keep_function : bool
            True to keep the function object, else return str
        Optional keyword input parameter is for internal use only
        and may prevent json serializability.
        """

        MagneticNetwork_dict = dict()
        # The class name is added to the dict for deserialisation purpose
        MagneticNetwork_dict["__class__"] = "MagneticNetwork"
        return MagneticNetwork_dict

    def copy(self):
        """Creates a deepcopy of the object"""

        # Handle deepcopy of all the properties
        # Creates new object of the same type with the copied properties
        obj_copy = type(self)()
        return obj_copy

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""
        pass
