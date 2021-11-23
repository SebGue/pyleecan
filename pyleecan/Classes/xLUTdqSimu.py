# -*- coding: utf-8 -*-
# File generated according to Generator/ClassesRef/Simulation/xLUTdqSimu.csv
# WARNING! All changes made in this file will be lost!
"""Method code available at https://github.com/Eomys/pyleecan/tree/master/pyleecan/Methods/Simulation/xLUTdqSimu
"""

from os import linesep
from sys import getsizeof
from logging import getLogger
from ._check import set_array, check_var, raise_
from ..Functions.get_logger import get_logger
from ..Functions.save import save
from ..Functions.copy import copy
from ..Functions.load import load_init_dict
from ..Functions.Load.import_class import import_class
from .Simu1 import Simu1

# Import all class method
# Try/catch to remove unnecessary dependencies in unused method
try:
    from ..Methods.Simulation.xLUTdqSimu.get_norm_machine import get_norm_machine
except ImportError as error:
    get_norm_machine = error

try:
    from ..Methods.Simulation.xLUTdqSimu.is_norm_machine import is_norm_machine
except ImportError as error:
    is_norm_machine = error

try:
    from ..Methods.Simulation.xLUTdqSimu.get_torque_norm import get_torque_norm
except ImportError as error:
    get_torque_norm = error

try:
    from ..Methods.Simulation.xLUTdqSimu.get_flux_norm import get_flux_norm
except ImportError as error:
    get_flux_norm = error

try:
    from ..Methods.Simulation.xLUTdqSimu.get_current_norm import get_current_norm
except ImportError as error:
    get_current_norm = error

try:
    from ..Methods.Simulation.xLUTdqSimu.get_VarLoadCurrent import get_VarLoadCurrent
except ImportError as error:
    get_VarLoadCurrent = error


from numpy import array, array_equal
from ._check import InitUnKnowClassError
from .Electrical import Electrical
from .Magnetics import Magnetics
from .Structural import Structural
from .Force import Force
from .Loss import Loss
from .Machine import Machine
from .Input import Input
from .VarSimu import VarSimu
from .Post import Post


class xLUTdqSimu(Simu1):
    """Object to simulate the flux and torque of a (normalized) machine"""

    VERSION = 1

    # Check ImportError to remove unnecessary dependencies in unused method
    # cf Methods.Simulation.xLUTdqSimu.get_norm_machine
    if isinstance(get_norm_machine, ImportError):
        get_norm_machine = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdqSimu method get_norm_machine: "
                    + str(get_norm_machine)
                )
            )
        )
    else:
        get_norm_machine = get_norm_machine
    # cf Methods.Simulation.xLUTdqSimu.is_norm_machine
    if isinstance(is_norm_machine, ImportError):
        is_norm_machine = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdqSimu method is_norm_machine: "
                    + str(is_norm_machine)
                )
            )
        )
    else:
        is_norm_machine = is_norm_machine
    # cf Methods.Simulation.xLUTdqSimu.get_torque_norm
    if isinstance(get_torque_norm, ImportError):
        get_torque_norm = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdqSimu method get_torque_norm: "
                    + str(get_torque_norm)
                )
            )
        )
    else:
        get_torque_norm = get_torque_norm
    # cf Methods.Simulation.xLUTdqSimu.get_flux_norm
    if isinstance(get_flux_norm, ImportError):
        get_flux_norm = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdqSimu method get_flux_norm: " + str(get_flux_norm)
                )
            )
        )
    else:
        get_flux_norm = get_flux_norm
    # cf Methods.Simulation.xLUTdqSimu.get_current_norm
    if isinstance(get_current_norm, ImportError):
        get_current_norm = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdqSimu method get_current_norm: "
                    + str(get_current_norm)
                )
            )
        )
    else:
        get_current_norm = get_current_norm
    # cf Methods.Simulation.xLUTdqSimu.get_VarLoadCurrent
    if isinstance(get_VarLoadCurrent, ImportError):
        get_VarLoadCurrent = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdqSimu method get_VarLoadCurrent: "
                    + str(get_VarLoadCurrent)
                )
            )
        )
    else:
        get_VarLoadCurrent = get_VarLoadCurrent
    # save and copy methods are available in all object
    save = save
    copy = copy
    # get_logger method is available in all object
    get_logger = get_logger

    def __init__(
        self,
        S_ref=None,
        fill_factor=0.4,
        S_max=None,
        S_d_max=None,
        S_q_max=None,
        Nr=None,
        Nrev=None,
        tM=20,
        simu_dict=None,
        elec=None,
        mag=None,
        struct=None,
        force=None,
        loss=None,
        name="",
        desc="",
        machine=-1,
        input=-1,
        logger_name="Pyleecan.Simulation",
        var_simu=None,
        postproc_list=-1,
        index=None,
        path_result=None,
        layer=None,
        layer_log_warn=None,
        init_dict=None,
        init_str=None,
    ):
        """Constructor of the class. Can be use in three ways :
        - __init__ (arg1 = 1, arg3 = 5) every parameters have name and default values
            for pyleecan type, -1 will call the default constructor
        - __init__ (init_dict = d) d must be a dictionary with property names as keys
        - __init__ (init_str = s) s must be a string
        s is the file path to load

        ndarray or list can be given for Vector and Matrix
        object or dict can be given for pyleecan Object"""

        if init_str is not None:  # Load from a file
            init_dict = load_init_dict(init_str)[1]
        if init_dict is not None:  # Initialisation by dict
            assert type(init_dict) is dict
            # Overwrite default value with init_dict content
            if "S_ref" in list(init_dict.keys()):
                S_ref = init_dict["S_ref"]
            if "fill_factor" in list(init_dict.keys()):
                fill_factor = init_dict["fill_factor"]
            if "S_max" in list(init_dict.keys()):
                S_max = init_dict["S_max"]
            if "S_d_max" in list(init_dict.keys()):
                S_d_max = init_dict["S_d_max"]
            if "S_q_max" in list(init_dict.keys()):
                S_q_max = init_dict["S_q_max"]
            if "Nr" in list(init_dict.keys()):
                Nr = init_dict["Nr"]
            if "Nrev" in list(init_dict.keys()):
                Nrev = init_dict["Nrev"]
            if "tM" in list(init_dict.keys()):
                tM = init_dict["tM"]
            if "simu_dict" in list(init_dict.keys()):
                simu_dict = init_dict["simu_dict"]
            if "elec" in list(init_dict.keys()):
                elec = init_dict["elec"]
            if "mag" in list(init_dict.keys()):
                mag = init_dict["mag"]
            if "struct" in list(init_dict.keys()):
                struct = init_dict["struct"]
            if "force" in list(init_dict.keys()):
                force = init_dict["force"]
            if "loss" in list(init_dict.keys()):
                loss = init_dict["loss"]
            if "name" in list(init_dict.keys()):
                name = init_dict["name"]
            if "desc" in list(init_dict.keys()):
                desc = init_dict["desc"]
            if "machine" in list(init_dict.keys()):
                machine = init_dict["machine"]
            if "input" in list(init_dict.keys()):
                input = init_dict["input"]
            if "logger_name" in list(init_dict.keys()):
                logger_name = init_dict["logger_name"]
            if "var_simu" in list(init_dict.keys()):
                var_simu = init_dict["var_simu"]
            if "postproc_list" in list(init_dict.keys()):
                postproc_list = init_dict["postproc_list"]
            if "index" in list(init_dict.keys()):
                index = init_dict["index"]
            if "path_result" in list(init_dict.keys()):
                path_result = init_dict["path_result"]
            if "layer" in list(init_dict.keys()):
                layer = init_dict["layer"]
            if "layer_log_warn" in list(init_dict.keys()):
                layer_log_warn = init_dict["layer_log_warn"]
        # Set the properties (value check and convertion are done in setter)
        self.S_ref = S_ref
        self.fill_factor = fill_factor
        self.S_max = S_max
        self.S_d_max = S_d_max
        self.S_q_max = S_q_max
        self.Nr = Nr
        self.Nrev = Nrev
        self.tM = tM
        self.simu_dict = simu_dict
        # Call Simu1 init
        super(xLUTdqSimu, self).__init__(
            elec=elec,
            mag=mag,
            struct=struct,
            force=force,
            loss=loss,
            name=name,
            desc=desc,
            machine=machine,
            input=input,
            logger_name=logger_name,
            var_simu=var_simu,
            postproc_list=postproc_list,
            index=index,
            path_result=path_result,
            layer=layer,
            layer_log_warn=layer_log_warn,
        )
        # The class is frozen (in Simu1 init), for now it's impossible to
        # add new properties

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        xLUTdqSimu_str = ""
        # Get the properties inherited from Simu1
        xLUTdqSimu_str += super(xLUTdqSimu, self).__str__()
        xLUTdqSimu_str += (
            "S_ref = "
            + linesep
            + str(self.S_ref).replace(linesep, linesep + "\t")
            + linesep
            + linesep
        )
        xLUTdqSimu_str += "fill_factor = " + str(self.fill_factor) + linesep
        xLUTdqSimu_str += "S_max = " + str(self.S_max) + linesep
        xLUTdqSimu_str += "S_d_max = " + str(self.S_d_max) + linesep
        xLUTdqSimu_str += "S_q_max = " + str(self.S_q_max) + linesep
        xLUTdqSimu_str += "Nr = " + str(self.Nr) + linesep
        xLUTdqSimu_str += "Nrev = " + str(self.Nrev) + linesep
        xLUTdqSimu_str += "tM = " + str(self.tM) + linesep
        xLUTdqSimu_str += "simu_dict = " + str(self.simu_dict) + linesep
        return xLUTdqSimu_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False

        # Check the properties inherited from Simu1
        if not super(xLUTdqSimu, self).__eq__(other):
            return False
        if not array_equal(other.S_ref, self.S_ref):
            return False
        if other.fill_factor != self.fill_factor:
            return False
        if other.S_max != self.S_max:
            return False
        if other.S_d_max != self.S_d_max:
            return False
        if other.S_q_max != self.S_q_max:
            return False
        if other.Nr != self.Nr:
            return False
        if other.Nrev != self.Nrev:
            return False
        if other.tM != self.tM:
            return False
        if other.simu_dict != self.simu_dict:
            return False
        return True

    def compare(self, other, name="self", ignore_list=None):
        """Compare two objects and return list of differences"""

        if ignore_list is None:
            ignore_list = list()
        if type(other) != type(self):
            return ["type(" + name + ")"]
        diff_list = list()

        # Check the properties inherited from Simu1
        diff_list.extend(super(xLUTdqSimu, self).compare(other, name=name))
        if not array_equal(other.S_ref, self.S_ref):
            diff_list.append(name + ".S_ref")
        if other._fill_factor != self._fill_factor:
            diff_list.append(name + ".fill_factor")
        if other._S_max != self._S_max:
            diff_list.append(name + ".S_max")
        if other._S_d_max != self._S_d_max:
            diff_list.append(name + ".S_d_max")
        if other._S_q_max != self._S_q_max:
            diff_list.append(name + ".S_q_max")
        if other._Nr != self._Nr:
            diff_list.append(name + ".Nr")
        if other._Nrev != self._Nrev:
            diff_list.append(name + ".Nrev")
        if other._tM != self._tM:
            diff_list.append(name + ".tM")
        if other._simu_dict != self._simu_dict:
            diff_list.append(name + ".simu_dict")
        # Filter ignore differences
        diff_list = list(filter(lambda x: x not in ignore_list, diff_list))
        return diff_list

    def __sizeof__(self):
        """Return the size in memory of the object (including all subobject)"""

        S = 0  # Full size of the object

        # Get size of the properties inherited from Simu1
        S += super(xLUTdqSimu, self).__sizeof__()
        S += getsizeof(self.S_ref)
        S += getsizeof(self.fill_factor)
        S += getsizeof(self.S_max)
        S += getsizeof(self.S_d_max)
        S += getsizeof(self.S_q_max)
        S += getsizeof(self.Nr)
        S += getsizeof(self.Nrev)
        S += getsizeof(self.tM)
        if self.simu_dict is not None:
            for key, value in self.simu_dict.items():
                S += getsizeof(value) + getsizeof(key)
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

        # Get the properties inherited from Simu1
        xLUTdqSimu_dict = super(xLUTdqSimu, self).as_dict(
            type_handle_ndarray=type_handle_ndarray,
            keep_function=keep_function,
            **kwargs
        )
        if self.S_ref is None:
            xLUTdqSimu_dict["S_ref"] = None
        else:
            if type_handle_ndarray == 0:
                xLUTdqSimu_dict["S_ref"] = self.S_ref.tolist()
            elif type_handle_ndarray == 1:
                xLUTdqSimu_dict["S_ref"] = self.S_ref.copy()
            elif type_handle_ndarray == 2:
                xLUTdqSimu_dict["S_ref"] = self.S_ref
            else:
                raise Exception(
                    "Unknown type_handle_ndarray: " + str(type_handle_ndarray)
                )
        xLUTdqSimu_dict["fill_factor"] = self.fill_factor
        xLUTdqSimu_dict["S_max"] = self.S_max
        xLUTdqSimu_dict["S_d_max"] = self.S_d_max
        xLUTdqSimu_dict["S_q_max"] = self.S_q_max
        xLUTdqSimu_dict["Nr"] = self.Nr
        xLUTdqSimu_dict["Nrev"] = self.Nrev
        xLUTdqSimu_dict["tM"] = self.tM
        xLUTdqSimu_dict["simu_dict"] = (
            self.simu_dict.copy() if self.simu_dict is not None else None
        )
        # The class name is added to the dict for deserialisation purpose
        # Overwrite the mother class name
        xLUTdqSimu_dict["__class__"] = "xLUTdqSimu"
        return xLUTdqSimu_dict

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""

        self.S_ref = None
        self.fill_factor = None
        self.S_max = None
        self.S_d_max = None
        self.S_q_max = None
        self.Nr = None
        self.Nrev = None
        self.tM = None
        self.simu_dict = None
        # Set to None the properties inherited from Simu1
        super(xLUTdqSimu, self)._set_None()

    def _get_S_ref(self):
        """getter of S_ref"""
        return self._S_ref

    def _set_S_ref(self, value):
        """setter of S_ref"""
        if type(value) is int and value == -1:
            value = array([])
        elif type(value) is list:
            try:
                value = array(value)
            except:
                pass
        check_var("S_ref", value, "ndarray")
        self._S_ref = value

    S_ref = property(
        fget=_get_S_ref,
        fset=_set_S_ref,
        doc=u"""Reference dq current density matrix

        :Type: ndarray
        """,
    )

    def _get_fill_factor(self):
        """getter of fill_factor"""
        return self._fill_factor

    def _set_fill_factor(self, value):
        """setter of fill_factor"""
        check_var("fill_factor", value, "float")
        self._fill_factor = value

    fill_factor = property(
        fget=_get_fill_factor,
        fset=_set_fill_factor,
        doc=u"""Reference winding fill factor

        :Type: float
        """,
    )

    def _get_S_max(self):
        """getter of S_max"""
        return self._S_max

    def _set_S_max(self, value):
        """setter of S_max"""
        check_var("S_max", value, "float")
        self._S_max = value

    S_max = property(
        fget=_get_S_max,
        fset=_set_S_max,
        doc=u"""Maximum rms current density

        :Type: float
        """,
    )

    def _get_S_d_max(self):
        """getter of S_d_max"""
        return self._S_d_max

    def _set_S_d_max(self, value):
        """setter of S_d_max"""
        check_var("S_d_max", value, "float")
        self._S_d_max = value

    S_d_max = property(
        fget=_get_S_d_max,
        fset=_set_S_d_max,
        doc=u"""Maximum d current density

        :Type: float
        """,
    )

    def _get_S_q_max(self):
        """getter of S_q_max"""
        return self._S_q_max

    def _set_S_q_max(self, value):
        """setter of S_q_max"""
        check_var("S_q_max", value, "float")
        self._S_q_max = value

    S_q_max = property(
        fget=_get_S_q_max,
        fset=_set_S_q_max,
        doc=u"""Maximum q current density

        :Type: float
        """,
    )

    def _get_Nr(self):
        """getter of Nr"""
        return self._Nr

    def _set_Nr(self, value):
        """setter of Nr"""
        check_var("Nr", value, "int", Vmin=1)
        self._Nr = value

    Nr = property(
        fget=_get_Nr,
        fset=_set_Nr,
        doc=u"""Rotation discretization

        :Type: int
        :min: 1
        """,
    )

    def _get_Nrev(self):
        """getter of Nrev"""
        return self._Nrev

    def _set_Nrev(self, value):
        """setter of Nrev"""
        check_var("Nrev", value, "float", Vmin=0)
        self._Nrev = value

    Nrev = property(
        fget=_get_Nrev,
        fset=_set_Nrev,
        doc=u"""Number of rotor revolution

        :Type: float
        :min: 0
        """,
    )

    def _get_tM(self):
        """getter of tM"""
        return self._tM

    def _set_tM(self, value):
        """setter of tM"""
        check_var("tM", value, "float")
        self._tM = value

    tM = property(
        fget=_get_tM,
        fset=_set_tM,
        doc=u"""Magnet temperature

        :Type: float
        """,
    )

    def _get_simu_dict(self):
        """getter of simu_dict"""
        return self._simu_dict

    def _set_simu_dict(self, value):
        """setter of simu_dict"""
        if type(value) is int and value == -1:
            value = dict()
        check_var("simu_dict", value, "dict")
        self._simu_dict = value

    simu_dict = property(
        fget=_get_simu_dict,
        fset=_set_simu_dict,
        doc=u"""Dict to set user MagFEMM parameter

        :Type: dict
        """,
    )
