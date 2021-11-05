# -*- coding: utf-8 -*-
# File generated according to Generator/ClassesRef/Simulation/MagLUT.csv
# WARNING! All changes made in this file will be lost!
"""Method code available at https://github.com/Eomys/pyleecan/tree/master/pyleecan/Methods/Simulation/MagLUT
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
from .XOutput import XOutput

# Import all class method
# Try/catch to remove unnecessary dependencies in unused method
try:
    from ..Methods.Simulation.MagLUT._set_fluxlinkage import _set_fluxlinkage
except ImportError as error:
    _set_fluxlinkage = error

try:
    from ..Methods.Simulation.MagLUT.get_interp import get_interp
except ImportError as error:
    get_interp = error

try:
    from ..Methods.Simulation.MagLUT.interpolate import interpolate
except ImportError as error:
    interpolate = error

try:
    from ..Methods.Simulation.MagLUT.comp_inductance import comp_inductance
except ImportError as error:
    comp_inductance = error

try:
    from ..Methods.Simulation.MagLUT.comp_grad_flux import comp_grad_flux
except ImportError as error:
    comp_grad_flux = error

try:
    from ..Methods.Simulation.MagLUT.get_fluxlinkage import get_fluxlinkage
except ImportError as error:
    get_fluxlinkage = error

try:
    from ..Methods.Simulation.MagLUT.get_torque import get_torque
except ImportError as error:
    get_torque = error


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


class MagLUT(XOutput):
    """Object to store the flux and torque of a (normalized) machine"""

    VERSION = 1

    # Check ImportError to remove unnecessary dependencies in unused method
    # cf Methods.Simulation.MagLUT._set_fluxlinkage
    if isinstance(_set_fluxlinkage, ImportError):
        _set_fluxlinkage = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagLUT method _set_fluxlinkage: " + str(_set_fluxlinkage)
                )
            )
        )
    else:
        _set_fluxlinkage = _set_fluxlinkage
    # cf Methods.Simulation.MagLUT.get_interp
    if isinstance(get_interp, ImportError):
        get_interp = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagLUT method get_interp: " + str(get_interp))
            )
        )
    else:
        get_interp = get_interp
    # cf Methods.Simulation.MagLUT.interpolate
    if isinstance(interpolate, ImportError):
        interpolate = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagLUT method interpolate: " + str(interpolate))
            )
        )
    else:
        interpolate = interpolate
    # cf Methods.Simulation.MagLUT.comp_inductance
    if isinstance(comp_inductance, ImportError):
        comp_inductance = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagLUT method comp_inductance: " + str(comp_inductance)
                )
            )
        )
    else:
        comp_inductance = comp_inductance
    # cf Methods.Simulation.MagLUT.comp_grad_flux
    if isinstance(comp_grad_flux, ImportError):
        comp_grad_flux = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagLUT method comp_grad_flux: " + str(comp_grad_flux)
                )
            )
        )
    else:
        comp_grad_flux = comp_grad_flux
    # cf Methods.Simulation.MagLUT.get_fluxlinkage
    if isinstance(get_fluxlinkage, ImportError):
        get_fluxlinkage = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagLUT method get_fluxlinkage: " + str(get_fluxlinkage)
                )
            )
        )
    else:
        get_fluxlinkage = get_fluxlinkage
    # cf Methods.Simulation.MagLUT.get_torque
    if isinstance(get_torque, ImportError):
        get_torque = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagLUT method get_torque: " + str(get_torque))
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
        # Call XOutput init
        super(MagLUT, self).__init__(
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
        # The class is frozen (in XOutput init), for now it's impossible to
        # add new properties

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        MagLUT_str = ""
        # Get the properties inherited from XOutput
        MagLUT_str += super(MagLUT, self).__str__()
        return MagLUT_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False

        # Check the properties inherited from XOutput
        if not super(MagLUT, self).__eq__(other):
            return False
        return True

    def compare(self, other, name="self", ignore_list=None):
        """Compare two objects and return list of differences"""

        if ignore_list is None:
            ignore_list = list()
        if type(other) != type(self):
            return ["type(" + name + ")"]
        diff_list = list()

        # Check the properties inherited from XOutput
        diff_list.extend(super(MagLUT, self).compare(other, name=name))
        # Filter ignore differences
        diff_list = list(filter(lambda x: x not in ignore_list, diff_list))
        return diff_list

    def __sizeof__(self):
        """Return the size in memory of the object (including all subobject)"""

        S = 0  # Full size of the object

        # Get size of the properties inherited from XOutput
        S += super(MagLUT, self).__sizeof__()
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

        # Get the properties inherited from XOutput
        MagLUT_dict = super(MagLUT, self).as_dict(
            type_handle_ndarray=type_handle_ndarray,
            keep_function=keep_function,
            **kwargs
        )
        # The class name is added to the dict for deserialisation purpose
        # Overwrite the mother class name
        MagLUT_dict["__class__"] = "MagLUT"
        return MagLUT_dict

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""

        # Set to None the properties inherited from XOutput
        super(MagLUT, self)._set_None()
