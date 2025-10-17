# -*- coding: utf-8 -*-
# File generated according to Generator/ClassesRef/Machine/MagnetWithGap.csv
# WARNING! All changes made in this file will be lost!
"""Method code available at https://github.com/Eomys/pyleecan/tree/master/pyleecan/Methods/Machine/MagnetWithGap
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
from .Magnet import Magnet

# Import all class method
# Try/catch to remove unnecessary dependencies in unused method
try:
    from ..Methods.Machine.MagnetWithGap.get_Hc import get_Hc
except ImportError as error:
    get_Hc = error

try:
    from ..Methods.Machine.MagnetWithGap.get_Brm import get_Brm
except ImportError as error:
    get_Brm = error

try:
    from ..Methods.Machine.MagnetWithGap.get_mur_lin import get_mur_lin
except ImportError as error:
    get_mur_lin = error

try:
    from ..Methods.Machine.MagnetWithGap.get_rho_elec import get_rho_elec
except ImportError as error:
    get_rho_elec = error


from numpy import isnan
from ._check import InitUnKnowClassError


class MagnetWithGap(Magnet):
    """Magnet for Hole Rotors with Gap"""

    VERSION = 1

    # Check ImportError to remove unnecessary dependencies in unused method
    # cf Methods.Machine.MagnetWithGap.get_Hc
    if isinstance(get_Hc, ImportError):
        get_Hc = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagnetWithGap method get_Hc: " + str(get_Hc))
            )
        )
    else:
        get_Hc = get_Hc
    # cf Methods.Machine.MagnetWithGap.get_Brm
    if isinstance(get_Brm, ImportError):
        get_Brm = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagnetWithGap method get_Brm: " + str(get_Brm))
            )
        )
    else:
        get_Brm = get_Brm
    # cf Methods.Machine.MagnetWithGap.get_mur_lin
    if isinstance(get_mur_lin, ImportError):
        get_mur_lin = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagnetWithGap method get_mur_lin: " + str(get_mur_lin)
                )
            )
        )
    else:
        get_mur_lin = get_mur_lin
    # cf Methods.Machine.MagnetWithGap.get_rho_elec
    if isinstance(get_rho_elec, ImportError):
        get_rho_elec = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagnetWithGap method get_rho_elec: " + str(get_rho_elec)
                )
            )
        )
    else:
        get_rho_elec = get_rho_elec
    # generic save method is available in all object
    save = save
    # get_logger method is available in all object
    get_logger = get_logger

    def __init__(
        self,
        hGap=0.0005,
        wGap=0.0005,
        mat_type=-1,
        type_magnetization=0,
        Lmag=0.95,
        Nseg=1,
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
            if "hGap" in list(init_dict.keys()):
                hGap = init_dict["hGap"]
            if "wGap" in list(init_dict.keys()):
                wGap = init_dict["wGap"]
            if "mat_type" in list(init_dict.keys()):
                mat_type = init_dict["mat_type"]
            if "type_magnetization" in list(init_dict.keys()):
                type_magnetization = init_dict["type_magnetization"]
            if "Lmag" in list(init_dict.keys()):
                Lmag = init_dict["Lmag"]
            if "Nseg" in list(init_dict.keys()):
                Nseg = init_dict["Nseg"]
        # Set the properties (value check and convertion are done in setter)
        self.hGap = hGap
        self.wGap = wGap
        # Call Magnet init
        super(MagnetWithGap, self).__init__(
            mat_type=mat_type,
            type_magnetization=type_magnetization,
            Lmag=Lmag,
            Nseg=Nseg,
        )
        # The class is frozen (in Magnet init), for now it's impossible to
        # add new properties

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        MagnetWithGap_str = ""
        # Get the properties inherited from Magnet
        MagnetWithGap_str += super(MagnetWithGap, self).__str__()
        MagnetWithGap_str += "hGap = " + str(self.hGap) + linesep
        MagnetWithGap_str += "wGap = " + str(self.wGap) + linesep
        return MagnetWithGap_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False

        # Check the properties inherited from Magnet
        if not super(MagnetWithGap, self).__eq__(other):
            return False
        if other.hGap != self.hGap:
            return False
        if other.wGap != self.wGap:
            return False
        return True

    def compare(self, other, name="self", ignore_list=None, is_add_value=False):
        """Compare two objects and return list of differences"""

        if ignore_list is None:
            ignore_list = list()
        if type(other) != type(self):
            return ["type(" + name + ")"]
        diff_list = list()

        # Check the properties inherited from Magnet
        diff_list.extend(
            super(MagnetWithGap, self).compare(
                other, name=name, ignore_list=ignore_list, is_add_value=is_add_value
            )
        )
        if (
            other._hGap is not None
            and self._hGap is not None
            and isnan(other._hGap)
            and isnan(self._hGap)
        ):
            pass
        elif other._hGap != self._hGap:
            if is_add_value:
                val_str = (
                    " (self=" + str(self._hGap) + ", other=" + str(other._hGap) + ")"
                )
                diff_list.append(name + ".hGap" + val_str)
            else:
                diff_list.append(name + ".hGap")
        if (
            other._wGap is not None
            and self._wGap is not None
            and isnan(other._wGap)
            and isnan(self._wGap)
        ):
            pass
        elif other._wGap != self._wGap:
            if is_add_value:
                val_str = (
                    " (self=" + str(self._wGap) + ", other=" + str(other._wGap) + ")"
                )
                diff_list.append(name + ".wGap" + val_str)
            else:
                diff_list.append(name + ".wGap")
        # Filter ignore differences
        diff_list = list(filter(lambda x: x not in ignore_list, diff_list))
        return diff_list

    def __sizeof__(self):
        """Return the size in memory of the object (including all subobject)"""

        S = 0  # Full size of the object

        # Get size of the properties inherited from Magnet
        S += super(MagnetWithGap, self).__sizeof__()
        S += getsizeof(self.hGap)
        S += getsizeof(self.wGap)
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

        # Get the properties inherited from Magnet
        MagnetWithGap_dict = super(MagnetWithGap, self).as_dict(
            type_handle_ndarray=type_handle_ndarray,
            keep_function=keep_function,
            **kwargs
        )
        MagnetWithGap_dict["hGap"] = self.hGap
        MagnetWithGap_dict["wGap"] = self.wGap
        # The class name is added to the dict for deserialisation purpose
        # Overwrite the mother class name
        MagnetWithGap_dict["__class__"] = "MagnetWithGap"
        return MagnetWithGap_dict

    def copy(self):
        """Creates a deepcopy of the object"""

        # Handle deepcopy of all the properties
        hGap_val = self.hGap
        wGap_val = self.wGap
        if self.mat_type is None:
            mat_type_val = None
        else:
            mat_type_val = self.mat_type.copy()
        type_magnetization_val = self.type_magnetization
        Lmag_val = self.Lmag
        Nseg_val = self.Nseg
        # Creates new object of the same type with the copied properties
        obj_copy = type(self)(
            hGap=hGap_val,
            wGap=wGap_val,
            mat_type=mat_type_val,
            type_magnetization=type_magnetization_val,
            Lmag=Lmag_val,
            Nseg=Nseg_val,
        )
        return obj_copy

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""

        self.hGap = None
        self.wGap = None
        # Set to None the properties inherited from Magnet
        super(MagnetWithGap, self)._set_None()

    def _get_hGap(self):
        """getter of hGap"""
        return self._hGap

    def _set_hGap(self, value):
        """setter of hGap"""
        check_var("hGap", value, "float", Vmin=0)
        self._hGap = value

    hGap = property(
        fget=_get_hGap,
        fset=_set_hGap,
        doc="""Gap to Lamination over Magnet Height

        :Type: float
        :min: 0
        """,
    )

    def _get_wGap(self):
        """getter of wGap"""
        return self._wGap

    def _set_wGap(self, value):
        """setter of wGap"""
        check_var("wGap", value, "float", Vmin=0)
        self._wGap = value

    wGap = property(
        fget=_get_wGap,
        fset=_set_wGap,
        doc="""Gap to Lamination over Magnet Width

        :Type: float
        :min: 0
        """,
    )
