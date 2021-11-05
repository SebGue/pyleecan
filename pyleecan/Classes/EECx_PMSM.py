# -*- coding: utf-8 -*-
# File generated according to Generator/ClassesRef/Simulation/EECx_PMSM.csv
# WARNING! All changes made in this file will be lost!
"""Method code available at https://github.com/Eomys/pyleecan/tree/master/pyleecan/Methods/Simulation/EECx_PMSM
"""

from os import linesep
from sys import getsizeof
from logging import getLogger
from ._check import check_var, raise_
from ..Functions.get_logger import get_logger
from ..Functions.save import save
from ..Functions.copy import copy
from ..Functions.load import load_init_dict
from ..Functions.Load.import_class import import_class
from .EECx import EECx

# Import all class method
# Try/catch to remove unnecessary dependencies in unused method
try:
    from ..Methods.Simulation.EECx_PMSM.solve import solve
except ImportError as error:
    solve = error

try:
    from ..Methods.Simulation.EECx_PMSM.comp_pmsm_model import comp_pmsm_model
except ImportError as error:
    comp_pmsm_model = error

try:
    from ..Methods.Simulation.EECx_PMSM.get_fun_core_loss import get_fun_core_loss
except ImportError as error:
    get_fun_core_loss = error

try:
    from ..Methods.Simulation.EECx_PMSM.get_fun_mech_loss import get_fun_mech_loss
except ImportError as error:
    get_fun_mech_loss = error

try:
    from ..Methods.Simulation.EECx_PMSM.get_torque_fun import get_torque_fun
except ImportError as error:
    get_torque_fun = error

try:
    from ..Methods.Simulation.EECx_PMSM.get_torque_fun import get_torque_fun
except ImportError as error:
    get_torque_fun = error

try:
    from ..Methods.Simulation.EECx_PMSM.get_fluxlinkage_fun import get_fluxlinkage_fun
except ImportError as error:
    get_fluxlinkage_fun = error

try:
    from ..Methods.Simulation.EECx_PMSM.comp_parameter import comp_parameter
except ImportError as error:
    comp_parameter = error

try:
    from ..Methods.Simulation.EECx_PMSM.get_fun_misc_loss import get_fun_misc_loss
except ImportError as error:
    get_fun_misc_loss = error


from ._check import InitUnKnowClassError
from .MagLUT import MagLUT
from .EndWindLeakage import EndWindLeakage
from .ACWinding import ACWinding
from .LossModel import LossModel
from .Machine import Machine


class EECx_PMSM(EECx):
    """Class of Electrical Equivalent Circuit for PMSM"""

    VERSION = 1

    # Check ImportError to remove unnecessary dependencies in unused method
    # cf Methods.Simulation.EECx_PMSM.solve
    if isinstance(solve, ImportError):
        solve = property(
            fget=lambda x: raise_(
                ImportError("Can't use EECx_PMSM method solve: " + str(solve))
            )
        )
    else:
        solve = solve
    # cf Methods.Simulation.EECx_PMSM.comp_pmsm_model
    if isinstance(comp_pmsm_model, ImportError):
        comp_pmsm_model = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use EECx_PMSM method comp_pmsm_model: "
                    + str(comp_pmsm_model)
                )
            )
        )
    else:
        comp_pmsm_model = comp_pmsm_model
    # cf Methods.Simulation.EECx_PMSM.get_fun_core_loss
    if isinstance(get_fun_core_loss, ImportError):
        get_fun_core_loss = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use EECx_PMSM method get_fun_core_loss: "
                    + str(get_fun_core_loss)
                )
            )
        )
    else:
        get_fun_core_loss = get_fun_core_loss
    # cf Methods.Simulation.EECx_PMSM.get_fun_mech_loss
    if isinstance(get_fun_mech_loss, ImportError):
        get_fun_mech_loss = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use EECx_PMSM method get_fun_mech_loss: "
                    + str(get_fun_mech_loss)
                )
            )
        )
    else:
        get_fun_mech_loss = get_fun_mech_loss
    # cf Methods.Simulation.EECx_PMSM.get_torque_fun
    if isinstance(get_torque_fun, ImportError):
        get_torque_fun = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use EECx_PMSM method get_torque_fun: " + str(get_torque_fun)
                )
            )
        )
    else:
        get_torque_fun = get_torque_fun
    # cf Methods.Simulation.EECx_PMSM.get_torque_fun
    if isinstance(get_torque_fun, ImportError):
        get_torque_fun = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use EECx_PMSM method get_torque_fun: " + str(get_torque_fun)
                )
            )
        )
    else:
        get_torque_fun = get_torque_fun
    # cf Methods.Simulation.EECx_PMSM.get_fluxlinkage_fun
    if isinstance(get_fluxlinkage_fun, ImportError):
        get_fluxlinkage_fun = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use EECx_PMSM method get_fluxlinkage_fun: "
                    + str(get_fluxlinkage_fun)
                )
            )
        )
    else:
        get_fluxlinkage_fun = get_fluxlinkage_fun
    # cf Methods.Simulation.EECx_PMSM.comp_parameter
    if isinstance(comp_parameter, ImportError):
        comp_parameter = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use EECx_PMSM method comp_parameter: " + str(comp_parameter)
                )
            )
        )
    else:
        comp_parameter = comp_parameter
    # cf Methods.Simulation.EECx_PMSM.get_fun_misc_loss
    if isinstance(get_fun_misc_loss, ImportError):
        get_fun_misc_loss = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use EECx_PMSM method get_fun_misc_loss: "
                    + str(get_fun_misc_loss)
                )
            )
        )
    else:
        get_fun_misc_loss = get_fun_misc_loss
    # save and copy methods are available in all object
    save = save
    copy = copy
    # get_logger method is available in all object
    get_logger = get_logger

    def __init__(
        self,
        table=-1,
        ew_leakage=None,
        winding=-1,
        core_loss=None,
        mech_loss=None,
        misc_loss=None,
        tW=20,
        machine=-1,
        logger_name="Pyleecan.Performance",
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
            if "table" in list(init_dict.keys()):
                table = init_dict["table"]
            if "ew_leakage" in list(init_dict.keys()):
                ew_leakage = init_dict["ew_leakage"]
            if "winding" in list(init_dict.keys()):
                winding = init_dict["winding"]
            if "core_loss" in list(init_dict.keys()):
                core_loss = init_dict["core_loss"]
            if "mech_loss" in list(init_dict.keys()):
                mech_loss = init_dict["mech_loss"]
            if "misc_loss" in list(init_dict.keys()):
                misc_loss = init_dict["misc_loss"]
            if "tW" in list(init_dict.keys()):
                tW = init_dict["tW"]
            if "machine" in list(init_dict.keys()):
                machine = init_dict["machine"]
            if "logger_name" in list(init_dict.keys()):
                logger_name = init_dict["logger_name"]
        # Set the properties (value check and convertion are done in setter)
        self.table = table
        self.ew_leakage = ew_leakage
        self.winding = winding
        self.core_loss = core_loss
        self.mech_loss = mech_loss
        self.misc_loss = misc_loss
        self.tW = tW
        # Call EECx init
        super(EECx_PMSM, self).__init__(machine=machine, logger_name=logger_name)
        # The class is frozen (in EECx init), for now it's impossible to
        # add new properties

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        EECx_PMSM_str = ""
        # Get the properties inherited from EECx
        EECx_PMSM_str += super(EECx_PMSM, self).__str__()
        if self.table is not None:
            tmp = self.table.__str__().replace(linesep, linesep + "\t").rstrip("\t")
            EECx_PMSM_str += "table = " + tmp
        else:
            EECx_PMSM_str += "table = None" + linesep + linesep
        if self.ew_leakage is not None:
            tmp = (
                self.ew_leakage.__str__().replace(linesep, linesep + "\t").rstrip("\t")
            )
            EECx_PMSM_str += "ew_leakage = " + tmp
        else:
            EECx_PMSM_str += "ew_leakage = None" + linesep + linesep
        if self.winding is not None:
            tmp = self.winding.__str__().replace(linesep, linesep + "\t").rstrip("\t")
            EECx_PMSM_str += "winding = " + tmp
        else:
            EECx_PMSM_str += "winding = None" + linesep + linesep
        if len(self.core_loss) == 0:
            EECx_PMSM_str += "core_loss = dict()" + linesep
        for key, obj in self.core_loss.items():
            tmp = (
                self.core_loss[key].__str__().replace(linesep, linesep + "\t") + linesep
            )
            EECx_PMSM_str += "core_loss[" + key + "] =" + tmp + linesep + linesep
        if len(self.mech_loss) == 0:
            EECx_PMSM_str += "mech_loss = dict()" + linesep
        for key, obj in self.mech_loss.items():
            tmp = (
                self.mech_loss[key].__str__().replace(linesep, linesep + "\t") + linesep
            )
            EECx_PMSM_str += "mech_loss[" + key + "] =" + tmp + linesep + linesep
        if len(self.misc_loss) == 0:
            EECx_PMSM_str += "misc_loss = dict()" + linesep
        for key, obj in self.misc_loss.items():
            tmp = (
                self.misc_loss[key].__str__().replace(linesep, linesep + "\t") + linesep
            )
            EECx_PMSM_str += "misc_loss[" + key + "] =" + tmp + linesep + linesep
        EECx_PMSM_str += "tW = " + str(self.tW) + linesep
        return EECx_PMSM_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False

        # Check the properties inherited from EECx
        if not super(EECx_PMSM, self).__eq__(other):
            return False
        if other.table != self.table:
            return False
        if other.ew_leakage != self.ew_leakage:
            return False
        if other.winding != self.winding:
            return False
        if other.core_loss != self.core_loss:
            return False
        if other.mech_loss != self.mech_loss:
            return False
        if other.misc_loss != self.misc_loss:
            return False
        if other.tW != self.tW:
            return False
        return True

    def compare(self, other, name="self", ignore_list=None):
        """Compare two objects and return list of differences"""

        if ignore_list is None:
            ignore_list = list()
        if type(other) != type(self):
            return ["type(" + name + ")"]
        diff_list = list()

        # Check the properties inherited from EECx
        diff_list.extend(super(EECx_PMSM, self).compare(other, name=name))
        if (other.table is None and self.table is not None) or (
            other.table is not None and self.table is None
        ):
            diff_list.append(name + ".table None mismatch")
        elif self.table is not None:
            diff_list.extend(self.table.compare(other.table, name=name + ".table"))
        if (other.ew_leakage is None and self.ew_leakage is not None) or (
            other.ew_leakage is not None and self.ew_leakage is None
        ):
            diff_list.append(name + ".ew_leakage None mismatch")
        elif self.ew_leakage is not None:
            diff_list.extend(
                self.ew_leakage.compare(other.ew_leakage, name=name + ".ew_leakage")
            )
        if (other.winding is None and self.winding is not None) or (
            other.winding is not None and self.winding is None
        ):
            diff_list.append(name + ".winding None mismatch")
        elif self.winding is not None:
            diff_list.extend(
                self.winding.compare(other.winding, name=name + ".winding")
            )
        if (other.core_loss is None and self.core_loss is not None) or (
            other.core_loss is not None and self.core_loss is None
        ):
            diff_list.append(name + ".core_loss None mismatch")
        elif self.core_loss is None:
            pass
        elif len(other.core_loss) != len(self.core_loss):
            diff_list.append("len(" + name + "core_loss)")
        else:
            for key in self.core_loss:
                diff_list.extend(
                    self.core_loss[key].compare(
                        other.core_loss[key], name=name + ".core_loss"
                    )
                )
        if (other.mech_loss is None and self.mech_loss is not None) or (
            other.mech_loss is not None and self.mech_loss is None
        ):
            diff_list.append(name + ".mech_loss None mismatch")
        elif self.mech_loss is None:
            pass
        elif len(other.mech_loss) != len(self.mech_loss):
            diff_list.append("len(" + name + "mech_loss)")
        else:
            for key in self.mech_loss:
                diff_list.extend(
                    self.mech_loss[key].compare(
                        other.mech_loss[key], name=name + ".mech_loss"
                    )
                )
        if (other.misc_loss is None and self.misc_loss is not None) or (
            other.misc_loss is not None and self.misc_loss is None
        ):
            diff_list.append(name + ".misc_loss None mismatch")
        elif self.misc_loss is None:
            pass
        elif len(other.misc_loss) != len(self.misc_loss):
            diff_list.append("len(" + name + "misc_loss)")
        else:
            for key in self.misc_loss:
                diff_list.extend(
                    self.misc_loss[key].compare(
                        other.misc_loss[key], name=name + ".misc_loss"
                    )
                )
        if other._tW != self._tW:
            diff_list.append(name + ".tW")
        # Filter ignore differences
        diff_list = list(filter(lambda x: x not in ignore_list, diff_list))
        return diff_list

    def __sizeof__(self):
        """Return the size in memory of the object (including all subobject)"""

        S = 0  # Full size of the object

        # Get size of the properties inherited from EECx
        S += super(EECx_PMSM, self).__sizeof__()
        S += getsizeof(self.table)
        S += getsizeof(self.ew_leakage)
        S += getsizeof(self.winding)
        if self.core_loss is not None:
            for key, value in self.core_loss.items():
                S += getsizeof(value) + getsizeof(key)
        if self.mech_loss is not None:
            for key, value in self.mech_loss.items():
                S += getsizeof(value) + getsizeof(key)
        if self.misc_loss is not None:
            for key, value in self.misc_loss.items():
                S += getsizeof(value) + getsizeof(key)
        S += getsizeof(self.tW)
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

        # Get the properties inherited from EECx
        EECx_PMSM_dict = super(EECx_PMSM, self).as_dict(
            type_handle_ndarray=type_handle_ndarray,
            keep_function=keep_function,
            **kwargs
        )
        if self.table is None:
            EECx_PMSM_dict["table"] = None
        else:
            EECx_PMSM_dict["table"] = self.table.as_dict(
                type_handle_ndarray=type_handle_ndarray,
                keep_function=keep_function,
                **kwargs
            )
        if self.ew_leakage is None:
            EECx_PMSM_dict["ew_leakage"] = None
        else:
            EECx_PMSM_dict["ew_leakage"] = self.ew_leakage.as_dict(
                type_handle_ndarray=type_handle_ndarray,
                keep_function=keep_function,
                **kwargs
            )
        if self.winding is None:
            EECx_PMSM_dict["winding"] = None
        else:
            EECx_PMSM_dict["winding"] = self.winding.as_dict(
                type_handle_ndarray=type_handle_ndarray,
                keep_function=keep_function,
                **kwargs
            )
        if self.core_loss is None:
            EECx_PMSM_dict["core_loss"] = None
        else:
            EECx_PMSM_dict["core_loss"] = dict()
            for key, obj in self.core_loss.items():
                if obj is not None:
                    EECx_PMSM_dict["core_loss"][key] = obj.as_dict(
                        type_handle_ndarray=type_handle_ndarray,
                        keep_function=keep_function,
                        **kwargs
                    )
                else:
                    EECx_PMSM_dict["core_loss"][key] = None
        if self.mech_loss is None:
            EECx_PMSM_dict["mech_loss"] = None
        else:
            EECx_PMSM_dict["mech_loss"] = dict()
            for key, obj in self.mech_loss.items():
                if obj is not None:
                    EECx_PMSM_dict["mech_loss"][key] = obj.as_dict(
                        type_handle_ndarray=type_handle_ndarray,
                        keep_function=keep_function,
                        **kwargs
                    )
                else:
                    EECx_PMSM_dict["mech_loss"][key] = None
        if self.misc_loss is None:
            EECx_PMSM_dict["misc_loss"] = None
        else:
            EECx_PMSM_dict["misc_loss"] = dict()
            for key, obj in self.misc_loss.items():
                if obj is not None:
                    EECx_PMSM_dict["misc_loss"][key] = obj.as_dict(
                        type_handle_ndarray=type_handle_ndarray,
                        keep_function=keep_function,
                        **kwargs
                    )
                else:
                    EECx_PMSM_dict["misc_loss"][key] = None
        EECx_PMSM_dict["tW"] = self.tW
        # The class name is added to the dict for deserialisation purpose
        # Overwrite the mother class name
        EECx_PMSM_dict["__class__"] = "EECx_PMSM"
        return EECx_PMSM_dict

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""

        if self.table is not None:
            self.table._set_None()
        if self.ew_leakage is not None:
            self.ew_leakage._set_None()
        if self.winding is not None:
            self.winding._set_None()
        self.core_loss = None
        self.mech_loss = None
        self.misc_loss = None
        self.tW = None
        # Set to None the properties inherited from EECx
        super(EECx_PMSM, self)._set_None()

    def _get_table(self):
        """getter of table"""
        return self._table

    def _set_table(self, value):
        """setter of table"""
        if isinstance(value, str):  # Load from file
            value = load_init_dict(value)[1]
        if isinstance(value, dict) and "__class__" in value:
            class_obj = import_class(
                "pyleecan.Classes", value.get("__class__"), "table"
            )
            value = class_obj(init_dict=value)
        elif type(value) is int and value == -1:  # Default constructor
            value = MagLUT()
        check_var("table", value, "MagLUT")
        self._table = value

        if self._table is not None:
            self._table.parent = self

    table = property(
        fget=_get_table,
        fset=_set_table,
        doc=u"""machine characteristics table

        :Type: MagLUT
        """,
    )

    def _get_ew_leakage(self):
        """getter of ew_leakage"""
        return self._ew_leakage

    def _set_ew_leakage(self, value):
        """setter of ew_leakage"""
        if isinstance(value, str):  # Load from file
            value = load_init_dict(value)[1]
        if isinstance(value, dict) and "__class__" in value:
            class_obj = import_class(
                "pyleecan.Classes", value.get("__class__"), "ew_leakage"
            )
            value = class_obj(init_dict=value)
        elif type(value) is int and value == -1:  # Default constructor
            value = EndWindLeakage()
        check_var("ew_leakage", value, "EndWindLeakage")
        self._ew_leakage = value

        if self._ew_leakage is not None:
            self._ew_leakage.parent = self

    ew_leakage = property(
        fget=_get_ew_leakage,
        fset=_set_ew_leakage,
        doc=u"""end winding leakage model

        :Type: EndWindLeakage
        """,
    )

    def _get_winding(self):
        """getter of winding"""
        return self._winding

    def _set_winding(self, value):
        """setter of winding"""
        if isinstance(value, str):  # Load from file
            value = load_init_dict(value)[1]
        if isinstance(value, dict) and "__class__" in value:
            class_obj = import_class(
                "pyleecan.Classes", value.get("__class__"), "winding"
            )
            value = class_obj(init_dict=value)
        elif type(value) is int and value == -1:  # Default constructor
            value = ACWinding()
        check_var("winding", value, "ACWinding")
        self._winding = value

        if self._winding is not None:
            self._winding.parent = self

    winding = property(
        fget=_get_winding,
        fset=_set_winding,
        doc=u"""model to compute the winding resistance

        :Type: ACWinding
        """,
    )

    def _get_core_loss(self):
        """getter of core_loss"""
        if self._core_loss is not None:
            for key, obj in self._core_loss.items():
                if obj is not None:
                    obj.parent = self
        return self._core_loss

    def _set_core_loss(self, value):
        """setter of core_loss"""
        if type(value) is dict:
            for key, obj in value.items():
                if type(obj) is dict:
                    class_obj = import_class(
                        "pyleecan.Classes", obj.get("__class__"), "core_loss"
                    )
                    value[key] = class_obj(init_dict=obj)
        if type(value) is int and value == -1:
            value = dict()
        check_var("core_loss", value, "{LossModel}")
        self._core_loss = value

    core_loss = property(
        fget=_get_core_loss,
        fset=_set_core_loss,
        doc=u"""dict of core loss models (e.g. magnets, frame and iron losses)

        :Type: {LossModel}
        """,
    )

    def _get_mech_loss(self):
        """getter of mech_loss"""
        if self._mech_loss is not None:
            for key, obj in self._mech_loss.items():
                if obj is not None:
                    obj.parent = self
        return self._mech_loss

    def _set_mech_loss(self, value):
        """setter of mech_loss"""
        if type(value) is dict:
            for key, obj in value.items():
                if type(obj) is dict:
                    class_obj = import_class(
                        "pyleecan.Classes", obj.get("__class__"), "mech_loss"
                    )
                    value[key] = class_obj(init_dict=obj)
        if type(value) is int and value == -1:
            value = dict()
        check_var("mech_loss", value, "{LossModel}")
        self._mech_loss = value

    mech_loss = property(
        fget=_get_mech_loss,
        fset=_set_mech_loss,
        doc=u"""dict of mechanical loss models (e.g. bearing and windage losses)

        :Type: {LossModel}
        """,
    )

    def _get_misc_loss(self):
        """getter of misc_loss"""
        if self._misc_loss is not None:
            for key, obj in self._misc_loss.items():
                if obj is not None:
                    obj.parent = self
        return self._misc_loss

    def _set_misc_loss(self, value):
        """setter of misc_loss"""
        if type(value) is dict:
            for key, obj in value.items():
                if type(obj) is dict:
                    class_obj = import_class(
                        "pyleecan.Classes", obj.get("__class__"), "misc_loss"
                    )
                    value[key] = class_obj(init_dict=obj)
        if type(value) is int and value == -1:
            value = dict()
        check_var("misc_loss", value, "{LossModel}")
        self._misc_loss = value

    misc_loss = property(
        fget=_get_misc_loss,
        fset=_set_misc_loss,
        doc=u"""dict of miscellaneous loss models

        :Type: {LossModel}
        """,
    )

    def _get_tW(self):
        """getter of tW"""
        return self._tW

    def _set_tW(self, value):
        """setter of tW"""
        check_var("tW", value, "float")
        self._tW = value

    tW = property(
        fget=_get_tW,
        fset=_set_tW,
        doc=u"""Stator winding temperature

        :Type: float
        """,
    )
