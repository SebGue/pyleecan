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
from .Magnetics import Magnetics

# Import all class method
# Try/catch to remove unnecessary dependencies in unused method
try:
    from ..Methods.Simulation.MagNetwork.assembler import assembler
except ImportError as error:
    assembler = error

try:
    from ..Methods.Simulation.MagNetwork.generalize_geometry import generalize_geometry
except ImportError as error:
    generalize_geometry = error

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
    from ..Methods.Simulation.MagNetwork.run_1 import run_1
except ImportError as error:
    run_1 = error

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

try:
    from ..Methods.Simulation.MagNetwork.comp_flux_airgap import comp_flux_airgap
except ImportError as error:
    comp_flux_airgap = error

try:
    from ..Methods.Simulation.MagNetwork.comp_flux_airgap_local import (
        comp_flux_airgap_local,
    )
except ImportError as error:
    comp_flux_airgap_local = error

try:
    from ..Methods.Simulation.MagNetwork.cartesianmeshclass_pyleecan import (
        cartesianmeshclass_pyleecan,
    )
except ImportError as error:
    cartesianmeshclass_pyleecan = error


from numpy import isnan
from ._check import InitUnKnowClassError


class MagNetwork(Magnetics):
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
    # cf Methods.Simulation.MagNetwork.generalize_geometry
    if isinstance(generalize_geometry, ImportError):
        generalize_geometry = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method generalize_geometry: "
                    + str(generalize_geometry)
                )
            )
        )
    else:
        generalize_geometry = generalize_geometry
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
    # cf Methods.Simulation.MagNetwork.run_1
    if isinstance(run_1, ImportError):
        run_1 = property(
            fget=lambda x: raise_(
                ImportError("Can't use MagNetwork method run_1: " + str(run_1))
            )
        )
    else:
        run_1 = run_1
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
    # cf Methods.Simulation.MagNetwork.comp_flux_airgap
    if isinstance(comp_flux_airgap, ImportError):
        comp_flux_airgap = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method comp_flux_airgap: "
                    + str(comp_flux_airgap)
                )
            )
        )
    else:
        comp_flux_airgap = comp_flux_airgap
    # cf Methods.Simulation.MagNetwork.comp_flux_airgap_local
    if isinstance(comp_flux_airgap_local, ImportError):
        comp_flux_airgap_local = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method comp_flux_airgap_local: "
                    + str(comp_flux_airgap_local)
                )
            )
        )
    else:
        comp_flux_airgap_local = comp_flux_airgap_local
    # cf Methods.Simulation.MagNetwork.cartesianmeshclass_pyleecan
    if isinstance(cartesianmeshclass_pyleecan, ImportError):
        cartesianmeshclass_pyleecan = property(
            fget=lambda x: raise_(
                ImportError(
                    "Can't use MagNetwork method cartesianmeshclass_pyleecan: "
                    + str(cartesianmeshclass_pyleecan)
                )
            )
        )
    else:
        cartesianmeshclass_pyleecan = cartesianmeshclass_pyleecan
    # generic save method is available in all object
    save = save
    # get_logger method is available in all object
    get_logger = get_logger

    def __init__(
        self,
        is_remove_slotS=False,
        is_remove_slotR=False,
        is_remove_ventS=False,
        is_remove_ventR=False,
        is_mmfs=True,
        is_mmfr=True,
        type_BH_stator=0,
        type_BH_rotor=0,
        is_periodicity_t=False,
        is_periodicity_a=False,
        angle_stator_shift=0,
        angle_rotor_shift=0,
        logger_name="Pyleecan.Magnetics",
        Slice_enforced=None,
        Nslices_enforced=None,
        type_distribution_enforced=None,
        is_current_harm=True,
        T_mag=20,
        is_periodicity_rotor=False,
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
            if "is_remove_slotS" in list(init_dict.keys()):
                is_remove_slotS = init_dict["is_remove_slotS"]
            if "is_remove_slotR" in list(init_dict.keys()):
                is_remove_slotR = init_dict["is_remove_slotR"]
            if "is_remove_ventS" in list(init_dict.keys()):
                is_remove_ventS = init_dict["is_remove_ventS"]
            if "is_remove_ventR" in list(init_dict.keys()):
                is_remove_ventR = init_dict["is_remove_ventR"]
            if "is_mmfs" in list(init_dict.keys()):
                is_mmfs = init_dict["is_mmfs"]
            if "is_mmfr" in list(init_dict.keys()):
                is_mmfr = init_dict["is_mmfr"]
            if "type_BH_stator" in list(init_dict.keys()):
                type_BH_stator = init_dict["type_BH_stator"]
            if "type_BH_rotor" in list(init_dict.keys()):
                type_BH_rotor = init_dict["type_BH_rotor"]
            if "is_periodicity_t" in list(init_dict.keys()):
                is_periodicity_t = init_dict["is_periodicity_t"]
            if "is_periodicity_a" in list(init_dict.keys()):
                is_periodicity_a = init_dict["is_periodicity_a"]
            if "angle_stator_shift" in list(init_dict.keys()):
                angle_stator_shift = init_dict["angle_stator_shift"]
            if "angle_rotor_shift" in list(init_dict.keys()):
                angle_rotor_shift = init_dict["angle_rotor_shift"]
            if "logger_name" in list(init_dict.keys()):
                logger_name = init_dict["logger_name"]
            if "Slice_enforced" in list(init_dict.keys()):
                Slice_enforced = init_dict["Slice_enforced"]
            if "Nslices_enforced" in list(init_dict.keys()):
                Nslices_enforced = init_dict["Nslices_enforced"]
            if "type_distribution_enforced" in list(init_dict.keys()):
                type_distribution_enforced = init_dict["type_distribution_enforced"]
            if "is_current_harm" in list(init_dict.keys()):
                is_current_harm = init_dict["is_current_harm"]
            if "T_mag" in list(init_dict.keys()):
                T_mag = init_dict["T_mag"]
            if "is_periodicity_rotor" in list(init_dict.keys()):
                is_periodicity_rotor = init_dict["is_periodicity_rotor"]
        # Set the properties (value check and convertion are done in setter)
        # Call Magnetics init
        super(MagNetwork, self).__init__(
            is_remove_slotS=is_remove_slotS,
            is_remove_slotR=is_remove_slotR,
            is_remove_ventS=is_remove_ventS,
            is_remove_ventR=is_remove_ventR,
            is_mmfs=is_mmfs,
            is_mmfr=is_mmfr,
            type_BH_stator=type_BH_stator,
            type_BH_rotor=type_BH_rotor,
            is_periodicity_t=is_periodicity_t,
            is_periodicity_a=is_periodicity_a,
            angle_stator_shift=angle_stator_shift,
            angle_rotor_shift=angle_rotor_shift,
            logger_name=logger_name,
            Slice_enforced=Slice_enforced,
            Nslices_enforced=Nslices_enforced,
            type_distribution_enforced=type_distribution_enforced,
            is_current_harm=is_current_harm,
            T_mag=T_mag,
            is_periodicity_rotor=is_periodicity_rotor,
        )
        # The class is frozen (in Magnetics init), for now it's impossible to
        # add new properties

    def __str__(self):
        """Convert this object in a readeable string (for print)"""

        MagNetwork_str = ""
        # Get the properties inherited from Magnetics
        MagNetwork_str += super(MagNetwork, self).__str__()
        return MagNetwork_str

    def __eq__(self, other):
        """Compare two objects (skip parent)"""

        if type(other) != type(self):
            return False

        # Check the properties inherited from Magnetics
        if not super(MagNetwork, self).__eq__(other):
            return False
        return True

    def compare(self, other, name="self", ignore_list=None, is_add_value=False):
        """Compare two objects and return list of differences"""

        if ignore_list is None:
            ignore_list = list()
        if type(other) != type(self):
            return ["type(" + name + ")"]
        diff_list = list()

        # Check the properties inherited from Magnetics
        diff_list.extend(
            super(MagNetwork, self).compare(
                other, name=name, ignore_list=ignore_list, is_add_value=is_add_value
            )
        )
        # Filter ignore differences
        diff_list = list(filter(lambda x: x not in ignore_list, diff_list))
        return diff_list

    def __sizeof__(self):
        """Return the size in memory of the object (including all subobject)"""

        S = 0  # Full size of the object

        # Get size of the properties inherited from Magnetics
        S += super(MagNetwork, self).__sizeof__()
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

        # Get the properties inherited from Magnetics
        MagNetwork_dict = super(MagNetwork, self).as_dict(
            type_handle_ndarray=type_handle_ndarray,
            keep_function=keep_function,
            **kwargs
        )
        # The class name is added to the dict for deserialisation purpose
        # Overwrite the mother class name
        MagNetwork_dict["__class__"] = "MagNetwork"
        return MagNetwork_dict

    def copy(self):
        """Creates a deepcopy of the object"""

        # Handle deepcopy of all the properties
        is_remove_slotS_val = self.is_remove_slotS
        is_remove_slotR_val = self.is_remove_slotR
        is_remove_ventS_val = self.is_remove_ventS
        is_remove_ventR_val = self.is_remove_ventR
        is_mmfs_val = self.is_mmfs
        is_mmfr_val = self.is_mmfr
        type_BH_stator_val = self.type_BH_stator
        type_BH_rotor_val = self.type_BH_rotor
        is_periodicity_t_val = self.is_periodicity_t
        is_periodicity_a_val = self.is_periodicity_a
        angle_stator_shift_val = self.angle_stator_shift
        angle_rotor_shift_val = self.angle_rotor_shift
        logger_name_val = self.logger_name
        if self.Slice_enforced is None:
            Slice_enforced_val = None
        else:
            Slice_enforced_val = self.Slice_enforced.copy()
        Nslices_enforced_val = self.Nslices_enforced
        type_distribution_enforced_val = self.type_distribution_enforced
        is_current_harm_val = self.is_current_harm
        T_mag_val = self.T_mag
        is_periodicity_rotor_val = self.is_periodicity_rotor
        # Creates new object of the same type with the copied properties
        obj_copy = type(self)(
            is_remove_slotS=is_remove_slotS_val,
            is_remove_slotR=is_remove_slotR_val,
            is_remove_ventS=is_remove_ventS_val,
            is_remove_ventR=is_remove_ventR_val,
            is_mmfs=is_mmfs_val,
            is_mmfr=is_mmfr_val,
            type_BH_stator=type_BH_stator_val,
            type_BH_rotor=type_BH_rotor_val,
            is_periodicity_t=is_periodicity_t_val,
            is_periodicity_a=is_periodicity_a_val,
            angle_stator_shift=angle_stator_shift_val,
            angle_rotor_shift=angle_rotor_shift_val,
            logger_name=logger_name_val,
            Slice_enforced=Slice_enforced_val,
            Nslices_enforced=Nslices_enforced_val,
            type_distribution_enforced=type_distribution_enforced_val,
            is_current_harm=is_current_harm_val,
            T_mag=T_mag_val,
            is_periodicity_rotor=is_periodicity_rotor_val,
        )
        return obj_copy

    def _set_None(self):
        """Set all the properties to None (except pyleecan object)"""

        # Set to None the properties inherited from Magnetics
        super(MagNetwork, self)._set_None()
