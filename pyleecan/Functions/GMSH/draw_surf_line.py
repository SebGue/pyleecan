from .get_boundary_condition import get_boundary_condition
from numpy import pi
from ...Classes.Arc import Arc
from ...Classes.Arc1 import Arc1
import cmath
from ...Functions.labels import BOUNDARY_PROP_LAB

tol = 1e-6


def draw_surf_line(
    surf,
    mesh_dict,
    boundary_prop,
    model,
    gmsh_dict,
    nsurf,
    mesh_size,
):
    """Draw the lines of a surface and handles the Arc>180deg

    Parameters
    ----------
    surf : Surface
        Surface object to draw
    mesh_dict : dict
        Dictionary to enforce the mesh (key: line id, value:nb element)
    boundary_prop : dict
        Dictionary to set the Boundary conditions
    model : Object
        Gmsh model
    gmsh_dict: dict
        dictionary containing the main parameters of GMSH File
    nsurf : int
        Index of the surface to draw
    mesh_size: float
        Default mesh element size

    Returns
    -------
    None
    """
    for ii, line in enumerate(surf.get_lines()):
        n_elem = None if not mesh_dict else mesh_dict[str(ii)]
        n_elem = n_elem if n_elem is not None else 0
        n_elem = n_elem if line.prop_dict else 0 # only element size if line has props
        bc_name = get_boundary_condition(line, boundary_prop)
        # Gmsh built-in engine does not allow arcs larger than 180deg
        # so arcs are split into two
        # TODO utilize OCC kernel to have 180° circles 
        if isinstance(line, Arc) and abs(line.get_angle() * 180.0 / pi) >= 180.0:
            rot_dir = 1 if line.is_trigo_direction == True else -1
            kwargs = dict(
                radius=rot_dir * line.comp_radius(),
                prop_dict=line.prop_dict,
                is_trigo_direction=line.is_trigo_direction,
            )
            arc1 = Arc1(begin=line.get_begin(), end=line.get_middle(), **kwargs)
            arc2 = Arc1(begin=line.get_middle(), end=line.get_end(), **kwargs)
            lines = [arc1, arc2]
        elif isinstance(line, Arc) and (abs(line.get_angle() * 180.0 / pi) <= tol):
            # Don't draw anything, this is a circle and usually is repeated ? TODO check
            lines = []
        # for debugging
        # if isinstance(line, Arc) and abs(line.get_angle() * 180.0 / pi) == 180.0: 
        #     lines = [line]
        else:
            lines = [line]

        for _line in lines:
            _add_line_to_dict(
                gmodel=model,
                line=_line,
                gmsh_dict=gmsh_dict,
                idx=nsurf,
                mesh_size=mesh_size,
                n_elements=n_elem,
                bc=bc_name,
            )


def _add_line_to_dict(
    gmodel, line, gmsh_dict={}, idx=0, mesh_size=1e-2, n_elements=0, bc=None
):
    """Draw a new line and add it to GMSH dictionary if it does not exist

    Parameters
    ----------
    gmodel : Object
        GMSH Model object
    line : Object
        Line Object
    d : Dictionary
        GMSH dictionary
    idx : int
        Surface index it belongs to
    mesh_size : float
        Points mesh size
    n_elements : int
        Number of elements on the line for meshing control
    bc : String
        Boundary condition name

    Returns
    -------
    None
    """

    ltag = None
    bz = line.get_begin()
    ez = line.get_end()
    bx, by = bz.real, bz.imag
    ex, ey = ez.real, ez.imag

    btag = gmodel.occ.addPoint(bx, by, 0, meshSize=mesh_size, tag=-1)
    etag = gmodel.occ.addPoint(ex, ey, 0, meshSize=mesh_size, tag=-1)

    if line.prop_dict and BOUNDARY_PROP_LAB in line.prop_dict:
        line_label = line.prop_dict[BOUNDARY_PROP_LAB]
    else:
        line_label = None
    if isinstance(line, Arc):
        cz = line.get_center()
        cx, cy = cz.real, cz.imag
        ctag = gmodel.occ.addPoint(cx, cy, 0, meshSize=mesh_size, tag=-1)
        ltag = gmodel.occ.addCircleArc(btag, ctag, etag, tag=-1)
        typ = "Circle"  # in terms of GMSH naming scheme
        arc_angle = cmath.phase(complex(ex, ey)) - cmath.phase(complex(bx, by))
        line_angle = None

    else:
        ltag = gmodel.occ.addLine(btag, etag, tag=-1)
        typ = "Line"
        ctag = None
        arc_angle = None
        line_angle = 0.5 * (cmath.phase(complex(ex, ey)) + cmath.phase(complex(bx, by)))

    gmodel.occ.synchronize()

    COM = gmodel.occ.getCenterOfMass(1, ltag)

    line_dict = {
        "tag": ltag,
        "label": line_label,
        "n_elements": n_elements,
        "bc_name": bc,
        "begin": {"tag": btag, "coord": (bx, by, 0)},
        "end": {"tag": etag, "coord": (ex, ey, 0)},
        "cent": {"tag": ctag, "coord": (cx, cy, 0)} if ctag else None,
        "arc_angle": arc_angle,
        "line_angle": line_angle,
        "COM": COM,
        "typ": typ,
    }

    gmsh_dict[idx]["lines"].append(line_dict)

    return None
