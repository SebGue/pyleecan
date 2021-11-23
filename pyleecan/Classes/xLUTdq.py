# -*- coding: utf-8 -*-
# File generated according to Generator/ClassesRef/Simulation/xLUTdq.csv
# WARNING! All changes made in this file will be lost!
"""Method code available at https://github.com/Eomys/pyleecan/tree/master/pyleecan/Methods/Simulation/xLUTdq
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
from .xLUT import xLUT

# Import all class method
# Try/catch to remove unnecessary dependencies in unused method
try:
    from ..Methods.Simulation.xLUTdq.get_param_dict import get_param_dict
except ImportError as error:
    get_param_dict = error

try:
    from ..Methods.Simulation.xLUTdq.get_bemf import get_bemf
except ImportError as error:
    get_bemf = error

try:
    from ..Methods.Simulation.xLUTdq.get_Ldqh import get_Ldqh
except ImportError as error:
    get_Ldqh = error

try:
    from ..Methods.Simulation.xLUTdq.get_Lmdqh import get_Lmdqh
except ImportError as error:
    get_Lmdqh = error

try:
    from ..Methods.Simulation.xLUTdq.import_from_data import import_from_data
except ImportError as error:
    import_from_data = error

try:
    from ..Methods.Simulation.xLUTdq.get_Phidqh_mean import get_Phidqh_mean
except ImportError as error:
    get_Phidqh_mean = error

try:
    from ..Methods.Simulation.xLUTdq.get_Phidqh_mag import get_Phidqh_mag
except ImportError as error:
    get_Phidqh_mag = error

try:
    from ..Methods.Simulation.xLUTdq.get_Phidqh_mag_mean import get_Phidqh_mag_mean
except ImportError as error:
    get_Phidqh_mag_mean = error

try:
    from ..Methods.Simulation.xLUTdq.get_Phidqh_mag_harm import get_Phidqh_mag_harm
except ImportError as error:
    get_Phidqh_mag_harm = error

try:
    from ..Methods.Simulation.xLUTdq.get_orders_dqh import get_orders_dqh
except ImportError as error:
    get_orders_dqh = error

try:
    from ..Methods.Simulation.xLUTdq.interp_Phi_dqh import interp_Phi_dqh
except ImportError as error:
    interp_Phi_dqh = error

try:
    from ..Methods.Simulation.xLUTdq._store_fluxlinkage import _store_fluxlinkage
except ImportError as error:
    _store_fluxlinkage = error

try:
    from ..Methods.Simulation.xLUTdq.get_interp import get_interp
except ImportError as error:
    get_interp = error

try:
    from ..Methods.Simulation.xLUTdq.interpolate import interpolate
except ImportError as error:
    interpolate = error

try:
    from ..Methods.Simulation.xLUTdq.comp_inductance import comp_inductance
except ImportError as error:
    comp_inductance = error

try:
    from ..Methods.Simulation.xLUTdq.comp_grad_flux import comp_grad_flux
except ImportError as error:
    comp_grad_flux = error

try:
    from ..Methods.Simulation.xLUTdq.get_fluxlinkage import get_fluxlinkage
except ImportError as error:
    get_fluxlinkage = error

try:
    from ..Methods.Simulation.xLUTdq.get_torque import get_torque
except ImportError as error:
    get_torque = error


from numpy import array, array_equal
from cloudpickle import dumps, loads
from ._check import CheckTypeError

try:
    from scipy.interpolate.interpolate import RegularGridInterpolator
except ImportError:
    RegularGridInterpolator = ImportError
from ._check import InitUnKnowClassError
from .ParamExplorer import ParamExplorer
from .Output import Output
from .DataKeeper import DataKeeper
from .Simulation import Simulation
from .OutGeo import OutGeo
from .OutElec import OutElec
from .OutMag import OutMag
from .OutStruct import OutStruct
from .OutPost import OutPost
from .OutForce import OutForce
from .OutLoss import OutLoss


class xLUTdq(xLUT):
    """Look Up Table class for dq OP matrix"""

    VERSION = 1

    # Check ImportError to remove unnecessary dependencies in unused method
    # cf Methods.Simulation.xLUTdq.get_param_dict
    if isinstance(get_param_dict, ImportError):
        get_param_dict = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method get_param_dict: " + str(get_param_dict)
                )
            )
        )
    else:
        get_param_dict = get_param_dict
    # cf Methods.Simulation.xLUTdq.get_bemf
    if isinstance(get_bemf, ImportError):
        get_bemf = property(
            fget=lambda x: raise_(
                ImportError("Can't use xLUTdq method get_bemf: " + str(get_bemf))
            )
        )
    else:
        get_bemf = get_bemf
    # cf Methods.Simulation.xLUTdq.get_Ldqh
    if isinstance(get_Ldqh, ImportError):
        get_Ldqh = property(
            fget=lambda x: raise_(
                ImportError("Can't use xLUTdq method get_Ldqh: " + str(get_Ldqh))
            )
        )
    else:
        get_Ldqh = get_Ldqh
    # cf Methods.Simulation.xLUTdq.get_Lmdqh
    if isinstance(get_Lmdqh, ImportError):
        get_Lmdqh = property(
            fget=lambda x: raise_(
                ImportError("Can't use xLUTdq method get_Lmdqh: " + str(get_Lmdqh))
            )
        )
    else:
        get_Lmdqh = get_Lmdqh
    # cf Methods.Simulation.xLUTdq.import_from_data
    if isinstance(import_from_data, ImportError):
        import_from_data = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method import_from_data: " + str(import_from_data)
                )
            )
        )
    else:
        import_from_data = import_from_data
    # cf Methods.Simulation.xLUTdq.get_Phidqh_mean
    if isinstance(get_Phidqh_mean, ImportError):
        get_Phidqh_mean = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method get_Phidqh_mean: " + str(get_Phidqh_mean)
                )
            )
        )
    else:
        get_Phidqh_mean = get_Phidqh_mean
    # cf Methods.Simulation.xLUTdq.get_Phidqh_mag
    if isinstance(get_Phidqh_mag, ImportError):
        get_Phidqh_mag = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method get_Phidqh_mag: " + str(get_Phidqh_mag)
                )
            )
        )
    else:
        get_Phidqh_mag = get_Phidqh_mag
    # cf Methods.Simulation.xLUTdq.get_Phidqh_mag_mean
    if isinstance(get_Phidqh_mag_mean, ImportError):
        get_Phidqh_mag_mean = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method get_Phidqh_mag_mean: "
                    + str(get_Phidqh_mag_mean)
                )
            )
        )
    else:
        get_Phidqh_mag_mean = get_Phidqh_mag_mean
    # cf Methods.Simulation.xLUTdq.get_Phidqh_mag_harm
    if isinstance(get_Phidqh_mag_harm, ImportError):
        get_Phidqh_mag_harm = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method get_Phidqh_mag_harm: "
                    + str(get_Phidqh_mag_harm)
                )
            )
        )
    else:
        get_Phidqh_mag_harm = get_Phidqh_mag_harm
    # cf Methods.Simulation.xLUTdq.get_orders_dqh
    if isinstance(get_orders_dqh, ImportError):
        get_orders_dqh = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method get_orders_dqh: " + str(get_orders_dqh)
                )
            )
        )
    else:
        get_orders_dqh = get_orders_dqh
    # cf Methods.Simulation.xLUTdq.interp_Phi_dqh
    if isinstance(interp_Phi_dqh, ImportError):
        interp_Phi_dqh = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method interp_Phi_dqh: " + str(interp_Phi_dqh)
                )
            )
        )
    else:
        interp_Phi_dqh = interp_Phi_dqh
    # cf Methods.Simulation.xLUTdq._store_fluxlinkage
    if isinstance(_store_fluxlinkage, ImportError):
        _store_fluxlinkage = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method _store_fluxlinkage: "
                    + str(_store_fluxlinkage)
                )
            )
        )
    else:
        _store_fluxlinkage = _store_fluxlinkage
    # cf Methods.Simulation.xLUTdq.get_interp
    if isinstance(get_interp, ImportError):
        get_interp = property(
            fget=lambda x: raise_(
                ImportError("Can't use xLUTdq method get_interp: " + str(get_interp))
            )
        )
    else:
        get_interp = get_interp
    # cf Methods.Simulation.xLUTdq.interpolate
    if isinstance(interpolate, ImportError):
        interpolate = property(
            fget=lambda x: raise_(
                ImportError("Can't use xLUTdq method interpolate: " + str(interpolate))
            )
        )
    else:
        interpolate = interpolate
    # cf Methods.Simulation.xLUTdq.comp_inductance
    if isinstance(comp_inductance, ImportError):
        comp_inductance = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method comp_inductance: " + str(comp_inductance)
                )
            )
        )
    else:
        comp_inductance = comp_inductance
    # cf Methods.Simulation.xLUTdq.comp_grad_flux
    if isinstance(comp_grad_flux, ImportError):
        comp_grad_flux = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method comp_grad_flux: " + str(comp_grad_flux)
                )
            )
        )
    else:
        comp_grad_flux = comp_grad_flux
    # cf Methods.Simulation.xLUTdq.get_fluxlinkage
    if isinstance(get_fluxlinkage, ImportError):
        get_fluxlinkage = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use xLUTdq method get_fluxlinkage: " + str(get_fluxlinkage)
                )
            )
        )
    else:
        get_fluxlinkage = get_fluxlinkage
    # cf Methods.Simulation.xLUTdq.get_torque
    if isinstance(get_torque, ImportError):
        get_torque = property(
            fget=lambda x: raise_(
                ImportError("Can't use xLUTdq method get_torque: " + str(get_torque))
            )
        )
    else:
        get_torque = get_torque
    # save and copy methods are available in all object
    save = save
    copy = copy
    # get_logger method is available in all object
    get_logger = get_logger

    def __init__(
        self,
        Phi_dqh_mean=None,
        Tmag_ref=20,
        Phi_dqh_mag=None,
        Phi_wind=None,
        Phi_dqh_interp=None,
        R1=None,
        L1=None,
        T1_ref=20,
        OP_matrix=None,
        phase_dir=None,
        paramexplorer_list=-1,
        output_list=-1,
        xoutput_dict=-1,
        nb_simu=0,
        xoutput_ref=None,
        xoutput_ref_index=None,
        simu=-1,
        path_result="",
        geo=-1,
        elec=-1,
        mag=-1,
        struct=-1,
        post=-1,
        logger_name="Pyleecan.Output",
        force=-1,
        loss=-1,
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
            if "Phi_dqh_mean" in list(init_dict.keys()):
                Phi_dqh_mean = init_dict["Phi_dqh_mean"]
            if "Tmag_ref" in list(init_dict.keys()):
                Tmag_ref = init_dict["Tmag_ref"]
            if "Phi_dqh_mag" in list(init_dict.keys()):
                Phi_dqh_mag = init_dict["Phi_dqh_mag"]
            if "Phi_wind" in list(init_dict.keys()):
                Phi_wind = init_dict["Phi_wind"]
            if "Phi_dqh_interp" in list(init_dict.keys()):
                Phi_dqh_interp = init_dict["Phi_dqh_interp"]
            if "R1" in list(init_dict.keys()):
                R1 = init_dict["R1"]
            if "L1" in list(init_dict.keys()):
                L1 = init_dict["L1"]
            if "T1_ref" in list(init_dict.keys()):
                T1_ref = init_dict["T1_ref"]
            if "OP_matrix" in list(init_dict.keys()):
                OP_matrix = init_dict["OP_matrix"]
            if "phase_dir" in list(init_dict.keys()):
                phase_dir = init_dict["phase_dir"]
            if "paramexplorer_list" in list(init_dict.keys()):
                paramexplorer_list = init_dict["paramexplorer_list"]
            if "output_list" in list(init_dict.keys()):
                output_list = init_dict["output_list"]
            if "xoutput_dict" in list(init_dict.keys()):
                xoutput_dict = init_dict["xoutput_dict"]
            if "nb_simu" in list(init_dict.keys()):
                nb_simu = init_dict["nb_simu"]
            if "xoutput_ref" in list(init_dict.keys()):
                xoutput_ref = init_dict["xoutput_ref"]
            if "xoutput_ref_index" in list(init_dict.keys()):
                xoutput_ref_index = init_dict["xoutput_ref_index"]
            if "simu" in list(init_dict.keys()):
                simu = init_dict["simu"]
            if "path_result" in list(init_dict.keys()):
                path_result = init_dict["path_result"]
            if "geo" in list(init_dict.keys()):
                geo = init_dict["geo"]
            if "elec" in list(init_dict.keys()):
                elec = init_dict["elec"]
            if "mag" in list(init_dict.keys()):
                mag = init_dict["mag"]
            if "struct" in list(init_dict.keys()):
                struct = init_dict["struct"]
            if "post" in list(init_dict.keys()):
                post = init_dict["post"]
            if "logger_name" in list(init_dict.keys()):
                logger_name = init_dict["logger_name"]
            if "force" in list(init_dict.keys()):
                force = init_dict["force"]
            if "loss" in list(init_dict.keys()):
                loss = init_dict["loss"]
        # Set the properties (value check and convertion are done in setter)
        self.Phi_dqh_mean = Phi_dqh_mean
        self.Tmag_ref = Tmag_ref
        self.Phi_dqh_mag = Phi_dqh_mag
        self.Phi_wind = Phi_wind
        self.Phi_dqh_interp = Phi_dqh_interp
        # Call xLUT init
        super(xLUTdq, self).__init__(
            R1=R1,
            L1=L1,
            T1_ref=T1_ref,
            OP_matrix=OP_matrix,
            phase_dir=phase_dir,
            paramexplorer_list=paramexplorer_list,
            output_list=output_list,
            xoutput_dict=xoutput_dict,
            nb_simu=nb_simu,
            xoutput_ref=xoutput_ref,
            xoutput_ref_index=xoutput_ref_index,
            simu=simu,
            path_result=path_result,
            geo=geo,
            elec=elec,
            mag=mag,
            struct=struct,
            post=post,
            logger_name=logger_name,
            force=force,
            loss=loss,
        )
        # The class is frozen (in xLUT init), for now it's impossible to
        # add new properties

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        xLUTdq_str = ""
        # Get the properties inherited from xLUT
        xLUTdq_str += super(xLUTdq, self).__str__()
        xLUTdq_str += (
            "Phi_dqh_mean = "
            + linesep
            + str(self.Phi_dqh_mean).replace(linesep, linesep + "\t")
            + linesep
            + linesep
        )
        xLUTdq_str += "Tmag_ref = " + str(self.Tmag_ref) + linesep
        xLUTdq_str += "Phi_dqh_mag = " + str(self.Phi_dqh_mag) + linesep + linesep
        xLUTdq_str += "Phi_wind = " + str(self.Phi_wind) + linesep + linesep
        xLUTdq_str += "Phi_dqh_interp = " + str(self.Phi_dqh_interp) + linesep + linesep
        return xLUTdq_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False

        # Check the properties inherited from xLUT
        if not super(xLUTdq, self).__eq__(other):
            return False
        if not array_equal(other.Phi_dqh_mean, self.Phi_dqh_mean):
            return False
        if other.Tmag_ref != self.Tmag_ref:
            return False
        if other.Phi_dqh_mag != self.Phi_dqh_mag:
            return False
        if other.Phi_wind != self.Phi_wind:
            return False
        if other.Phi_dqh_interp != self.Phi_dqh_interp:
            return False
        return True

    def compare(self, other, name="self", ignore_list=None):
        """Compare two objects and return list of differences"""

        if ignore_list is None:
            ignore_list = list()
        if type(other) != type(self):
            return ["type(" + name + ")"]
        diff_list = list()

        # Check the properties inherited from xLUT
        diff_list.extend(super(xLUTdq, self).compare(other, name=name))
        if not array_equal(other.Phi_dqh_mean, self.Phi_dqh_mean):
            diff_list.append(name + ".Phi_dqh_mean")
        if other._Tmag_ref != self._Tmag_ref:
            diff_list.append(name + ".Tmag_ref")
        if (other.Phi_dqh_mag is None and self.Phi_dqh_mag is not None) or (
            other.Phi_dqh_mag is not None and self.Phi_dqh_mag is None
        ):
            diff_list.append(name + ".Phi_dqh_mag None mismatch")
        elif self.Phi_dqh_mag is not None:
            diff_list.extend(
                self.Phi_dqh_mag.compare(other.Phi_dqh_mag, name=name + ".Phi_dqh_mag")
            )
        if (other.Phi_wind is None and self.Phi_wind is not None) or (
            other.Phi_wind is not None and self.Phi_wind is None
        ):
            diff_list.append(name + ".Phi_wind None mismatch")
        elif self.Phi_wind is None:
            pass
        elif len(other.Phi_wind) != len(self.Phi_wind):
            diff_list.append("len(" + name + ".Phi_wind)")
        else:
            for ii in range(len(other.Phi_wind)):
                diff_list.extend(
                    self.Phi_wind[ii].compare(
                        other.Phi_wind[ii], name=name + ".Phi_wind[" + str(ii) + "]"
                    )
                )
        if (other.Phi_dqh_interp is None and self.Phi_dqh_interp is not None) or (
            other.Phi_dqh_interp is not None and self.Phi_dqh_interp is None
        ):
            diff_list.append(name + ".Phi_dqh_interp None mismatch")
        elif (
            self.Phi_dqh_interp is not None
            and self.Phi_dqh_interp != other.Phi_dqh_interp
        ):
            diff_list.append(name + ".Phi_dqh_interp")
        # Filter ignore differences
        diff_list = list(filter(lambda x: x not in ignore_list, diff_list))
        return diff_list

    def __sizeof__(self):
        """Return the size in memory of the object (including all subobject)"""

        S = 0  # Full size of the object

        # Get size of the properties inherited from xLUT
        S += super(xLUTdq, self).__sizeof__()
        S += getsizeof(self.Phi_dqh_mean)
        S += getsizeof(self.Tmag_ref)
        S += getsizeof(self.Phi_dqh_mag)
        if self.Phi_wind is not None:
            for value in self.Phi_wind:
                S += getsizeof(value)
        S += getsizeof(self.Phi_dqh_interp)
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

        # Get the properties inherited from xLUT
        xLUTdq_dict = super(xLUTdq, self).as_dict(
            type_handle_ndarray=type_handle_ndarray,
            keep_function=keep_function,
            **kwargs
        )
        if self.Phi_dqh_mean is None:
            xLUTdq_dict["Phi_dqh_mean"] = None
        else:
            if type_handle_ndarray == 0:
                xLUTdq_dict["Phi_dqh_mean"] = self.Phi_dqh_mean.tolist()
            elif type_handle_ndarray == 1:
                xLUTdq_dict["Phi_dqh_mean"] = self.Phi_dqh_mean.copy()
            elif type_handle_ndarray == 2:
                xLUTdq_dict["Phi_dqh_mean"] = self.Phi_dqh_mean
            else:
                raise Exception(
                    "Unknown type_handle_ndarray: " + str(type_handle_ndarray)
                )
        xLUTdq_dict["Tmag_ref"] = self.Tmag_ref
        if self.Phi_dqh_mag is None:
            xLUTdq_dict["Phi_dqh_mag"] = None
        else:
            xLUTdq_dict["Phi_dqh_mag"] = self.Phi_dqh_mag.as_dict(
                type_handle_ndarray=type_handle_ndarray,
                keep_function=keep_function,
                **kwargs
            )
        if self.Phi_wind is None:
            xLUTdq_dict["Phi_wind"] = None
        else:
            xLUTdq_dict["Phi_wind"] = list()
            for obj in self.Phi_wind:
                if obj is not None:
                    xLUTdq_dict["Phi_wind"].append(
                        obj.as_dict(
                            type_handle_ndarray=type_handle_ndarray,
                            keep_function=keep_function,
                            **kwargs
                        )
                    )
                else:
                    xLUTdq_dict["Phi_wind"].append(None)
        if self.Phi_dqh_interp is None:
            xLUTdq_dict["Phi_dqh_interp"] = None
        else:
            # Store serialized data (using cloudpickle) and str
            # to read it in json save files
            xLUTdq_dict["Phi_dqh_interp"] = {
                "__class__": str(type(self._Phi_dqh_interp)),
                "__repr__": str(self._Phi_dqh_interp.__repr__()),
                "serialized": dumps(self._Phi_dqh_interp).decode("ISO-8859-2"),
            }
        # The class name is added to the dict for deserialisation purpose
        # Overwrite the mother class name
        xLUTdq_dict["__class__"] = "xLUTdq"
        return xLUTdq_dict

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""

        self.Phi_dqh_mean = None
        self.Tmag_ref = None
        self.Phi_dqh_mag = None
        self.Phi_wind = None
        self.Phi_dqh_interp = None
        # Set to None the properties inherited from xLUT
        super(xLUTdq, self)._set_None()

    def _get_Phi_dqh_mean(self):
        """getter of Phi_dqh_mean"""
        return self._Phi_dqh_mean

    def _set_Phi_dqh_mean(self, value):
        """setter of Phi_dqh_mean"""
        if type(value) is int and value == -1:
            value = array([])
        elif type(value) is list:
            try:
                value = array(value)
            except:
                pass
        check_var("Phi_dqh_mean", value, "ndarray")
        self._Phi_dqh_mean = value

    Phi_dqh_mean = property(
        fget=_get_Phi_dqh_mean,
        fset=_set_Phi_dqh_mean,
        doc=u"""RMS stator winding flux table in dqh frame (including magnets and currents given by I_dqh)

        :Type: ndarray
        """,
    )

    def _get_Tmag_ref(self):
        """getter of Tmag_ref"""
        return self._Tmag_ref

    def _set_Tmag_ref(self, value):
        """setter of Tmag_ref"""
        check_var("Tmag_ref", value, "float")
        self._Tmag_ref = value

    Tmag_ref = property(
        fget=_get_Tmag_ref,
        fset=_set_Tmag_ref,
        doc=u"""Magnet average temperature at which Phi_dqh is given

        :Type: float
        """,
    )

    def _get_Phi_dqh_mag(self):
        """getter of Phi_dqh_mag"""
        return self._Phi_dqh_mag

    def _set_Phi_dqh_mag(self, value):
        """setter of Phi_dqh_mag"""
        if isinstance(value, str):  # Load from file
            value = load_init_dict(value)[1]
        if isinstance(value, dict) and "__class__" in value:
            class_obj = import_class(
                "SciDataTool.Classes", value.get("__class__"), "Phi_dqh_mag"
            )
            value = class_obj(init_dict=value)
        elif type(value) is int and value == -1:  # Default constructor
            value = DataND()
        check_var("Phi_dqh_mag", value, "DataND")
        self._Phi_dqh_mag = value

    Phi_dqh_mag = property(
        fget=_get_Phi_dqh_mag,
        fset=_set_Phi_dqh_mag,
        doc=u"""RMS stator winding flux linkage spectrum in dqh frame including harmonics (only magnets)

        :Type: SciDataTool.Classes.DataND.DataND
        """,
    )

    def _get_Phi_wind(self):
        """getter of Phi_wind"""
        if self._Phi_wind is not None:
            for obj in self._Phi_wind:
                if obj is not None:
                    obj.parent = self
        return self._Phi_wind

    def _set_Phi_wind(self, value):
        """setter of Phi_wind"""
        if type(value) is list:
            for ii, obj in enumerate(value):
                if type(obj) is dict:
                    class_obj = import_class(
                        "SciDataTool.Classes", obj.get("__class__"), "Phi_wind"
                    )
                    value[ii] = class_obj(init_dict=obj)
                if value[ii] is not None:
                    value[ii].parent = self
        if value == -1:
            value = list()
        check_var("Phi_wind", value, "[DataND]")
        self._Phi_wind = value

    Phi_wind = property(
        fget=_get_Phi_wind,
        fset=_set_Phi_wind,
        doc=u"""Stator winding flux function of time and phases

        :Type: [SciDataTool.Classes.DataND.DataND]
        """,
    )

    def _get_Phi_dqh_interp(self):
        """getter of Phi_dqh_interp"""
        return self._Phi_dqh_interp

    def _set_Phi_dqh_interp(self, value):
        """setter of Phi_dqh_interp"""
        check_var("Phi_dqh_interp", value, "RegularGridInterpolator")
        self._Phi_dqh_interp = value

    Phi_dqh_interp = property(
        fget=_get_Phi_dqh_interp,
        fset=_set_Phi_dqh_interp,
        doc=u"""Interpolant function of Phi_dqh

        :Type: scipy.interpolate.interpolate.RegularGridInterpolator
        """,
    )
