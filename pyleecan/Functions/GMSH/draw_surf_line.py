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
        bc_name = get_boundary_condition(line, boundary_prop)
        # Gmsh built-in engine does not allow arcs larger than 180deg
        # so arcs are split into two
        if isinstance(line, Arc) and abs(line.get_angle() * 180.0 / pi) >= 180.0:
            rot_dir = 1 if line.is_trigo_direction == True else -1
            arc1 = Arc1(
                begin=line.get_begin(),
                end=line.get_middle(),
                radius=rot_dir * line.comp_radius(),
                prop_dict=line.prop_dict,
                is_trigo_direction=line.is_trigo_direction,
            )
            arc2 = Arc1(
                begin=line.get_middle(),
                end=line.get_end(),
                radius=rot_dir * line.comp_radius(),
                prop_dict=line.prop_dict,
                is_trigo_direction=line.is_trigo_direction,
            )
            for arc in [arc1, arc2]:
                _add_line_to_dict(
                    gmodel=model,
                    line=arc,
                    gmsh_dict=gmsh_dict,
                    idx=nsurf,
                    mesh_size=mesh_size,
                    n_elements=n_elem,
                    bc=bc_name,
                )
        elif isinstance(line, Arc) and (abs(line.get_angle() * 180.0 / pi) <= tol):
            # Don't draw anything, this is a circle and usually is repeated ?
            pass
        else:
            _add_line_to_dict(
                gmodel=model,
                line=line,
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

    dlines = list()
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
        if len(dlines) > 0:
            for iline in dlines:
                p = _find_points_from_line(d, iline)
                if p[0] == btag and p[1] == etag and p[2] == ctag:
                    ltag = iline
                    break
                elif p[0] == etag and p[1] == btag and p[2] == ctag:
                    ltag = -iline
                    break
                else:
                    pass
            if ltag is None:
                ltag = gmodel.occ.addCircleArc(btag, ctag, etag, tag=-1)
                if n_elements > 0:
                    gmodel.occ.synchronize()
                    gmodel.mesh.setTransfiniteCurve(ltag, n_elements + 1, "Progression")
        else:
            ltag = gmodel.occ.addCircleArc(btag, ctag, etag, tag=-1)
            if n_elements > 0:
                gmodel.occ.synchronize()
                gmodel.mesh.setTransfiniteCurve(ltag, n_elements + 1, "Progression")

        # To avoid fill the dictionary with repeated lines
        repeated = False
        # for lvalues in gmsh_dict[idx].values():
        #     if type(lvalues) is not dict:
        #         continue
        #     else:
        #         if lvalues["tag"] == ltag:
        #             repeated = True

        if not repeated:
            arc_angle = cmath.phase(complex(ex, ey)) - cmath.phase(complex(bx, by))
            gmsh_dict[idx]["lines"].append(
                {
                    "tag": ltag,
                    "label": line_label,
                    "n_elements": n_elements,
                    "bc_name": bc,
                    "begin": {"tag": btag, "coord": complex(bx, by)},
                    "end": {"tag": etag, "coord": complex(ex, ey)},
                    "cent": {"tag": ctag, "coord": complex(cx, cy)},
                    "arc_angle": arc_angle,
                    "line_angle": None,
                }
            )

    else:
        if len(dlines) > 0:
            for iline in dlines:
                p = _find_points_from_line(d, iline)
                if p[0] == btag and p[1] == etag:
                    ltag = iline
                    break
                elif p[0] == etag and p[1] == btag:
                    ltag = -iline
                    break
                else:
                    pass
            if ltag is None:
                ltag = gmodel.occ.addLine(btag, etag, tag=-1)
                if n_elements > 0:
                    gmodel.occ.synchronize()
                    gmodel.mesh.setTransfiniteCurve(ltag, n_elements + 1, "Progression")
        else:
            ltag = gmodel.occ.addLine(btag, etag, tag=-1)
            if n_elements > 0:
                gmodel.occ.synchronize()
                gmodel.mesh.setTransfiniteCurve(ltag, n_elements + 1, "Progression")

        # To avoid fill the dictionary with repeated lines
        repeated = False
        # for lvalues in gmsh_dict[idx].values():
        #     if type(lvalues) is not dict:
        #         continue
        #     else:
        #         if lvalues["tag"] == ltag:
        #             repeated = True

        if not repeated:
            line_angle = 0.5 * (
                cmath.phase(complex(ex, ey)) + cmath.phase(complex(bx, by))
            )
            gmsh_dict[idx]["lines"].append(
                {
                    "tag": ltag,
                    "label": line_label,
                    "n_elements": n_elements,
                    "bc_name": bc,
                    "begin": {"tag": btag, "coord": complex(bx, by)},
                    "end": {"tag": etag, "coord": complex(ex, ey)},
                    "arc_angle": None,
                    "line_angle": line_angle,
                }
            )

    return None


def _find_points_from_line(gmsh_dict={}, ltag=-1):
    """Find points tag from existing lines

    Parameters
    ----------
    d : Dictionary
        GMSH dictionary
    ltag : int
        line tag

    Returns
    -------
    coord : float
        Coordinates of point if found
    """
    btag = None
    etag = None
    ctag = None
    for surf in gmsh_dict.values():
        for lvalues in surf.values():
            if type(lvalues) is not dict:
                continue
            if lvalues["tag"] != ltag:
                continue
            else:
                for pid, pvalues in lvalues.items():
                    if type(pvalues) is not dict:
                        continue
                    if pid == "begin":
                        btag = pvalues["tag"]
                    elif pid == "end":
                        etag = pvalues["tag"]
                    elif pid == "cent":
                        ctag = pvalues["tag"]
                break
    return [btag, etag, ctag]
