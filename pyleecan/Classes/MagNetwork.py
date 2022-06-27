# -*- coding: utf-8 -*-
# File generated according to Generator/ClassesRef/Simulation/MagNetwork.csv
# WARNING! All changes made in this file will be lost!
"""Method code available at https://github.com/Eomys/pyleecan/tree/master/pyleecan/Methods/Simulation/MagNetwork
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
    from ..Methods.Simulation.MagNetwork.assembler import assembler
except ImportError as error:
    assembler = error

try:
    from ..Methods.Simulation.MagNetwork.geometry_linear_motor import (
        geometry_linear_motor,
    )
except ImportError as error:
    geometry_linear_motor = error

try:
    from ..Methods.Simulation.MagNetwork.geometry_linear_motor_separetion import (
        geometry_linear_motor_separetion,
    )
except ImportError as error:
    geometry_linear_motor_separetion = error

try:
    from ..Methods.Simulation.MagNetwork.plot.Plot_B import Plot_B
except ImportError as error:
    Plot_B = error

try:
    from ..Methods.Simulation.MagNetwork.plot.view_contour_flux import view_contour_flux
except ImportError as error:
    view_contour_flux = error

try:
    from ..Methods.Simulation.MagNetwork.plot.view_contour_flux2 import (
        view_contour_flux2,
    )
except ImportError as error:
    view_contour_flux2 = error

try:
    from ..Methods.Simulation.MagNetwork.post_processing import post_processing
except ImportError as error:
    post_processing = error

try:
    from ..Methods.Simulation.MagNetwork.run import run
except ImportError as error:
    run = error

try:
    from ..Methods.Simulation.MagNetwork.run_non_linear import run_non_linear
except ImportError as error:
    run_non_linear = error

try:
    from ..Methods.Simulation.MagNetwork.run_radial import run_radial
except ImportError as error:
    run_radial = error

try:
    from ..Methods.Simulation.MagNetwork.solver_linear_model import solver_linear_model
except ImportError as error:
    solver_linear_model = error

try:
    from ..Methods.Simulation.MagNetwork.solver_non_linear_model import (
        solver_non_linear_model,
    )
except ImportError as error:
    solver_non_linear_model = error

try:
    from ..Methods.Simulation.MagNetwork.pre_processing.init_cell import init_cell
except ImportError as error:
    init_cell = error

try:
    from ..Methods.Simulation.MagNetwork.pre_processing.init_mesh_BC import init_mesh_BC
except ImportError as error:
    init_mesh_BC = error

try:
    from ..Methods.Simulation.MagNetwork.pre_processing.init_permeabilty_cell import (
        init_permeabilty_cell,
    )
except ImportError as error:
    init_permeabilty_cell = error

try:
    from ..Methods.Simulation.MagNetwork.pre_processing.init_point import init_point
except ImportError as error:
    init_point = error

try:
    from ..Methods.Simulation.MagNetwork.pre_processing.init_reluc import init_reluc
except ImportError as error:
    init_reluc = error

try:
    from ..Methods.Simulation.MagNetwork.pre_processing.numeroting_unknows import (
        numeroting_unknows,
    )
except ImportError as error:
    numeroting_unknows = error

try:
    from ..Methods.Simulation.MagNetwork.pre_processing.save_mesh import save_mesh
except ImportError as error:
    save_mesh = error

try:
    from ..Methods.Simulation.MagNetwork.assembler.assembly import assembly
except ImportError as error:
    assembly = error

try:
    from ..Methods.Simulation.MagNetwork.assembler.assembly_one_area import (
        assembly_one_area,
    )
except ImportError as error:
    assembly_one_area = error

try:
    from ..Methods.Simulation.MagNetwork.assembler.right_member_assembly import (
        right_member_assembly,
    )
except ImportError as error:
    right_member_assembly = error

try:
    from ..Methods.Simulation.MagNetwork.post_processing.add_BC_to_F import add_BC_to_F
except ImportError as error:
    add_BC_to_F = error

try:
    from ..Methods.Simulation.MagNetwork.post_processing.compute_B_radial import (
        compute_B_radial,
    )
except ImportError as error:
    compute_B_radial = error

try:
    from ..Methods.Simulation.MagNetwork.post_processing.compute_B_square import (
        compute_B_square,
    )
except ImportError as error:
    compute_B_square = error


from numpy import isnan
from ._check import InitUnKnowClassError


class MagNetwork(FrozenClass):
    """Abstract class to solve the electric machines using the reluctance network"""

    VERSION = 1

    # Check ImportError to remove unnecessary dependencies in unused method
    # cf Methods.Simulation.MagNetwork.assembler
    if isinstance(assembler, ImportError):
        assembler = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagNetwork method assembler: " + str(assembler))
            )
        )
    else:
        assembler = assembler
    # cf Methods.Simulation.MagNetwork.geometry_linear_motor
    if isinstance(geometry_linear_motor, ImportError):
        geometry_linear_motor = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method geometry_linear_motor: "
                    + str(geometry_linear_motor)
                )
            )
        )
    else:
        geometry_linear_motor = geometry_linear_motor
    # cf Methods.Simulation.MagNetwork.geometry_linear_motor_separetion
    if isinstance(geometry_linear_motor_separetion, ImportError):
        geometry_linear_motor_separetion = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method geometry_linear_motor_separetion: "
                    + str(geometry_linear_motor_separetion)
                )
            )
        )
    else:
        geometry_linear_motor_separetion = geometry_linear_motor_separetion
    # cf Methods.Simulation.MagNetwork.plot.Plot_B
    if isinstance(Plot_B, ImportError):
        Plot_B = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagNetwork method Plot_B: " + str(Plot_B))
            )
        )
    else:
        Plot_B = Plot_B
    # cf Methods.Simulation.MagNetwork.plot.view_contour_flux
    if isinstance(view_contour_flux, ImportError):
        view_contour_flux = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method view_contour_flux: "
                    + str(view_contour_flux)
                )
            )
        )
    else:
        view_contour_flux = view_contour_flux
    # cf Methods.Simulation.MagNetwork.plot.view_contour_flux2
    if isinstance(view_contour_flux2, ImportError):
        view_contour_flux2 = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method view_contour_flux2: "
                    + str(view_contour_flux2)
                )
            )
        )
    else:
        view_contour_flux2 = view_contour_flux2
    # cf Methods.Simulation.MagNetwork.post_processing
    if isinstance(post_processing, ImportError):
        post_processing = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method post_processing: "
                    + str(post_processing)
                )
            )
        )
    else:
        post_processing = post_processing
    # cf Methods.Simulation.MagNetwork.run
    if isinstance(run, ImportError):
        run = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagNetwork method run: " + str(run))
            )
        )
    else:
        run = run
    # cf Methods.Simulation.MagNetwork.run_non_linear
    if isinstance(run_non_linear, ImportError):
        run_non_linear = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method run_non_linear: " + str(run_non_linear)
                )
            )
        )
    else:
        run_non_linear = run_non_linear
    # cf Methods.Simulation.MagNetwork.run_radial
    if isinstance(run_radial, ImportError):
        run_radial = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method run_radial: " + str(run_radial)
                )
            )
        )
    else:
        run_radial = run_radial
    # cf Methods.Simulation.MagNetwork.solver_linear_model
    if isinstance(solver_linear_model, ImportError):
        solver_linear_model = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method solver_linear_model: "
                    + str(solver_linear_model)
                )
            )
        )
    else:
        solver_linear_model = solver_linear_model
    # cf Methods.Simulation.MagNetwork.solver_non_linear_model
    if isinstance(solver_non_linear_model, ImportError):
        solver_non_linear_model = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method solver_non_linear_model: "
                    + str(solver_non_linear_model)
                )
            )
        )
    else:
        solver_non_linear_model = solver_non_linear_model
    # cf Methods.Simulation.MagNetwork.pre_processing.init_cell
    if isinstance(init_cell, ImportError):
        init_cell = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagNetwork method init_cell: " + str(init_cell))
            )
        )
    else:
        init_cell = init_cell
    # cf Methods.Simulation.MagNetwork.pre_processing.init_mesh_BC
    if isinstance(init_mesh_BC, ImportError):
        init_mesh_BC = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method init_mesh_BC: " + str(init_mesh_BC)
                )
            )
        )
    else:
        init_mesh_BC = init_mesh_BC
    # cf Methods.Simulation.MagNetwork.pre_processing.init_permeabilty_cell
    if isinstance(init_permeabilty_cell, ImportError):
        init_permeabilty_cell = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method init_permeabilty_cell: "
                    + str(init_permeabilty_cell)
                )
            )
        )
    else:
        init_permeabilty_cell = init_permeabilty_cell
    # cf Methods.Simulation.MagNetwork.pre_processing.init_point
    if isinstance(init_point, ImportError):
        init_point = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method init_point: " + str(init_point)
                )
            )
        )
    else:
        init_point = init_point
    # cf Methods.Simulation.MagNetwork.pre_processing.init_reluc
    if isinstance(init_reluc, ImportError):
        init_reluc = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method init_reluc: " + str(init_reluc)
                )
            )
        )
    else:
        init_reluc = init_reluc
    # cf Methods.Simulation.MagNetwork.pre_processing.numeroting_unknows
    if isinstance(numeroting_unknows, ImportError):
        numeroting_unknows = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method numeroting_unknows: "
                    + str(numeroting_unknows)
                )
            )
        )
    else:
        numeroting_unknows = numeroting_unknows
    # cf Methods.Simulation.MagNetwork.pre_processing.save_mesh
    if isinstance(save_mesh, ImportError):
        save_mesh = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagNetwork method save_mesh: " + str(save_mesh))
            )
        )
    else:
        save_mesh = save_mesh
    # cf Methods.Simulation.MagNetwork.assembler.assembly
    if isinstance(assembly, ImportError):
        assembly = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagNetwork method assembly: " + str(assembly))
            )
        )
    else:
        assembly = assembly
    # cf Methods.Simulation.MagNetwork.assembler.assembly_one_area
    if isinstance(assembly_one_area, ImportError):
        assembly_one_area = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method assembly_one_area: "
                    + str(assembly_one_area)
                )
            )
        )
    else:
        assembly_one_area = assembly_one_area
    # cf Methods.Simulation.MagNetwork.assembler.right_member_assembly
    if isinstance(right_member_assembly, ImportError):
        right_member_assembly = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method right_member_assembly: "
                    + str(right_member_assembly)
                )
            )
        )
    else:
        right_member_assembly = right_member_assembly
    # cf Methods.Simulation.MagNetwork.post_processing.add_BC_to_F
    if isinstance(add_BC_to_F, ImportError):
        add_BC_to_F = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method add_BC_to_F: " + str(add_BC_to_F)
                )
            )
        )
    else:
        add_BC_to_F = add_BC_to_F
    # cf Methods.Simulation.MagNetwork.post_processing.compute_B_radial
    if isinstance(compute_B_radial, ImportError):
        compute_B_radial = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method compute_B_radial: "
                    + str(compute_B_radial)
                )
            )
        )
    else:
        compute_B_radial = compute_B_radial
    # cf Methods.Simulation.MagNetwork.post_processing.compute_B_square
    if isinstance(compute_B_square, ImportError):
        compute_B_square = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method compute_B_square: "
                    + str(compute_B_square)
                )
            )
        )
    else:
        compute_B_square = compute_B_square
    # generic save method is available in all object
    save = save
    # get_logger method is available in all object
    get_logger = get_logger

    def __init__(self, file_path=-1, init_dict=None, init_str=None):
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
            if "file_path" in list(init_dict.keys()):
                file_path = init_dict["file_path"]
        # Set the properties (value check and convertion are done in setter)
        self.parent = None
        self.file_path = file_path

        # The class is frozen, for now it's impossible to add new properties
        self._freeze()

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        MagNetwork_str = ""
        if self.parent is None:
            MagNetwork_str += "parent = None " + linesep
        else:
            MagNetwork_str += "parent = " + str(type(self.parent)) + " object" + linesep
        MagNetwork_str += "file_path = " + str(self.file_path) + linesep + linesep
        return MagNetwork_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False
        if isinstance(self.file_path, np.ndarray) and not np.array_equal(
            other.file_path, self.file_path
        ):
            return False
        elif other.file_path != self.file_path:
            return False
        return True

    def compare(self, other, name="self", ignore_list=None, is_add_value=False):
        """Compare two objects and return list of differences"""

        if ignore_list is None:
            ignore_list = list()
        if type(other) != type(self):
            return ["type(" + name + ")"]
        diff_list = list()
        if (other.file_path is None and self.file_path is not None) or (
            other.file_path is not None and self.file_path is None
        ):
            diff_list.append(name + ".file_path")
        elif self.file_path is None:
            pass
        elif isinstance(self.file_path, np.ndarray) and not np.array_equal(
            other.file_path, self.file_path
        ):
            diff_list.append(name + ".file_path")
        elif hasattr(self.file_path, "compare"):
            diff_list.extend(
                self.file_path.compare(
                    other.file_path,
                    name=name + ".file_path",
                    ignore_list=ignore_list,
                    is_add_value=is_add_value,
                )
            )
        elif other._file_path != self._file_path:
            diff_list.append(name + ".file_path")
        # Filter ignore differences
        diff_list = list(filter(lambda x: x not in ignore_list, diff_list))
        return diff_list

    def __sizeof__(self):
        """Return the size in memory of the object (including all subobject)"""

        S = 0  # Full size of the object
        S += getsizeof(self.file_path)
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

        MagNetwork_dict = dict()
        if self.file_path is None:
            MagNetwork_dict["file_path"] = None
        elif isinstance(self.file_path, np.ndarray):
            if type_handle_ndarray == 0:
                MagNetwork_dict["file_path"] = self.file_path.tolist()
            elif type_handle_ndarray == 1:
                MagNetwork_dict["file_path"] = self.file_path.copy()
            elif type_handle_ndarray == 2:
                MagNetwork_dict["file_path"] = self.file_path
            else:
                raise Exception(
                    "Unknown type_handle_ndarray: " + str(type_handle_ndarray)
                )
        elif hasattr(self.file_path, "as_dict"):
            MagNetwork_dict["file_path"] = self.file_path.as_dict(
                type_handle_ndarray=type_handle_ndarray,
                keep_function=keep_function,
                **kwargs
            )
        else:
            MagNetwork_dict["file_path"] = self.file_path
        # The class name is added to the dict for deserialisation purpose
        MagNetwork_dict["__class__"] = "MagNetwork"
        return MagNetwork_dict

    def copy(self):
        """Creates a deepcopy of the object"""

        # Handle deepcopy of all the properties
        if hasattr(self.file_path, "copy"):
            file_path_val = self.file_path.copy()
        else:
            file_path_val = self.file_path
        # Creates new object of the same type with the copied properties
        obj_copy = type(self)(file_path=file_path_val)
        return obj_copy

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""

        if hasattr(self.file_path, "_set_None"):
            self.file_path._set_None()
        else:
            self.file_path = None

    def _get_file_path(self):
        """getter of file_path"""
        return self._file_path

    def _set_file_path(self, value):
        """setter of file_path"""
        if isinstance(value, dict) and "__class__" in value:
            try:
                class_obj = import_class(
                    "pyleecan.Classes", value.get("__class__"), "file_path"
                )
            except:
                class_obj = import_class(
                    "SciDataTool.Classes", value.get("__class__"), "file_path"
                )
            value = class_obj(init_dict=value)
        elif type(value) is list:
            try:
                value = np.array(value)
            except:
                pass
        check_var("file_path", value, "")
        self._file_path = value

        if hasattr(self._file_path, "parent"):
            self._file_path.parent = self

    file_path = property(
        fget=_get_file_path,
        fset=_set_file_path,
        doc=u"""file path of the .json file of the defined machine

        :Type: 
        """,
    )
