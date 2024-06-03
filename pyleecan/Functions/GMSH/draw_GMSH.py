from ...Functions.labels import (
    short_label,
    LAM_LAB_S,
    ROTOR_LAB_S,
    STATOR_LAB_S,
    decode_label,
    SLID_LAB,
    AIRGAP_LAB,
    BOT_LAB,
    TOP_LAB,
    AIRBOX_LAB,
    AIRBOX_R_LAB,
    AR_B_LAB,
    AR_T_LAB,
    SBR_B_LAB,
    SBR_T_LAB,
    AS_BL_LAB,
    AS_BR_LAB,
    AS_TL_LAB,
    AS_TR_LAB,
)
from ...Functions.GMSH import InputError
from ...Classes.SurfRing import SurfRing
from ...Functions.GMSH.get_sliding_band import get_sliding_band
from ...Functions.GMSH.get_air_box import get_air_box
from ...Functions.GMSH.draw_surf_line import draw_surf_line
from ...Functions.GMSH.comp_gmsh_mesh_dict import comp_gmsh_mesh_dict
import gmsh

from os import replace
from os.path import splitext

from numpy import pi, angle

from ...Functions.get_logger import get_logger


def draw_GMSH(
    output,
    sym,
    boundary_prop,
    is_lam_only_S=False,
    is_lam_only_R=False,
    user_mesh_dict={},
    path_save="GMSH_model.msh",
    is_sliding_band=False,
    is_airbox=False,
    is_set_labels=False,
    is_run=False,
):
    """Draws a machine mesh in GMSH format

    Parameters
    ----------
    output : Output
        Output object
    sym : int
        the symmetry applied on the stator and the rotor (take into account antiperiodicity)
    boundary_prop : dict
        dictionary to match FEA boundary conditions (dict values) with line boundary property (dict keys)
    is_lam_only_S: bool
        Draw only stator lamination
    is_lam_only_R: bool
        Draw only rotor lamination
    user_mesh_dict :dict
        Dictionary to enforce the mesh size on some surface/lines (key: surface label, value dict:key=line id, value:nb element)
    path_save : str
        Path to save the result msh file
    is_sliding_band : bool
        True uses sliding band, else airgap
    is_airbox : bool
        True to add the airbox
    is_set_label : bool
        True to set all line labels as physical groups
    is_run : bool
        True to launch Gmsh GUI at the end

    Returns
    -------
    GMSH_dict : dict
        dictionary containing the main parameters of GMSH File
    """
    # check some input parameter
    if is_lam_only_S and is_lam_only_R:
        raise InputError(
            "is_lam_only_S and is_lam_only_R can't be True at the same time"
        )

    # get machine
    machine = output.simu.machine
    mesh_dict = {}

    # Default stator mesh element size
    mesh_size_S = machine.stator.Rext / 100.0  # Stator
    mesh_size_R = machine.rotor.Rext / 25.0  # Rotor
    mesh_size_SB = 2.0 * pi * machine.rotor.Rext / 360.0  # Sliding Band
    mesh_size_AB = machine.stator.Rext / 50.0  # AirBox
    lam_list = machine.get_lam_list()
    lam_int = lam_list[0]
    lam_ext = lam_list[1]
    lab_int = lam_int.get_label()
    lab_ext = lam_ext.get_label()

    # For readibility
    model = gmsh.model
    factory = model.occ

    # Start a new model
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", int(False))
    gmsh.option.setNumber("Geometry.CopyMeshingMethod", 1)
    gmsh.option.setNumber("Geometry.PointNumbers", 0)
    gmsh.option.setNumber("Geometry.LineNumbers", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min(mesh_size_S, mesh_size_R))
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", max(mesh_size_S, mesh_size_R))
    model.add("Pyleecan")

    alpha = 0

    #####################
    # Adding Rotor
    #####################
    gmsh_dict = {}
    nSurf = 0  # number of surfaces
    # Drawing Rotor and Shaft surfaces
    if not is_lam_only_S:
        # get Pyleecan rotor surface definitions
        surf_list = []
        if machine.shaft is not None:
            surf_list.extend(machine.shaft.build_geometry(sym=sym, alpha=alpha))
        surf_list.extend(machine.rotor.build_geometry(sym=sym, alpha=alpha))

        # draw all surfaces lines in gmsh and store needed information in gmsh_dict
        for surf in surf_list:
            nSurf += 1  # surface number
            args = [boundary_prop, mesh_size_R, user_mesh_dict]
            _draw_lines(model, surf, gmsh_dict, nSurf, *args)

        # draw surfaces (from lines) in gmsh and add surface tags to gmsh_dict
        for surf in gmsh_dict.values():
            _draw_surf(factory, surf)

        # get lamination index
        idLam = [id for id, s in gmsh_dict.items() if LAM_LAB_S in s["label"]][0]

        # cut all holes from lamination
        tools = [(2, surf["tag"]) for idx, surf in gmsh_dict.items() if idx != idLam]
        surf = gmsh_dict[idLam]
        surf["tag"] = model.occ.cut(
            [(2, surf["tag"])], tools, removeObject=True, removeTool=False
        )[0][0][1]

        factory.synchronize()

    # store rotor dict
    rotor_dict = gmsh_dict.copy()

    #####################
    # Adding Stator
    #####################
    gmsh_dict = {}  # init new dict for stator

    # nSurf = 0
    if not is_lam_only_R:
        stator_list = list()
        stator_list.extend(machine.stator.build_geometry(sym=sym, alpha=alpha))

        for surf in stator_list:
            nSurf += 1
            args = [boundary_prop, mesh_size_S, user_mesh_dict]
            _draw_lines(model, surf, gmsh_dict, nSurf, *args)

        # draw surfaces (from lines) in gmsh and add surface tags to gmsh_dict
        for surf in gmsh_dict.values():
            _draw_surf(factory, surf)

        # get lamination index
        idLam = [id for id, s in gmsh_dict.items() if LAM_LAB_S in s["label"]][0]

        # cut all holes from lamination
        tools = [(2, surf["tag"]) for idx, surf in gmsh_dict.items() if idx != idLam]
        surf = gmsh_dict[idLam]
        surf["tag"] = model.occ.cut(
            [(2, surf["tag"])], tools, removeObject=True, removeTool=False
        )[0][0][1]

        factory.synchronize()

    stator_dict = gmsh_dict.copy()

    #####################
    # Adding Sliding Band
    #####################
    gmsh_dict = {}

    if is_sliding_band and (not is_lam_only_R) and (not is_lam_only_S):
        sb_list = get_sliding_band(sym=sym, machine=machine)

        for surf in sb_list:
            nSurf += 1
            args = [boundary_prop, mesh_size_SB, user_mesh_dict]
            _draw_lines(model, surf, gmsh_dict, nSurf, *args)

        # draw surfaces (from lines) in gmsh and add surface tags to gmsh_dict
        for surf in gmsh_dict.values():
            _draw_surf(factory, surf)

        # get indizes of lamination air gaps
        for idx, surf in gmsh_dict.items():
            label = surf["label"]
            if ROTOR_LAB_S in label and AIRGAP_LAB in label:
                idRotorAG = idx  # index of lamination
            if STATOR_LAB_S in label and AIRGAP_LAB in label:
                idStatorAG = idx  # index of lamination

        # cut airgap regions
        rotor_tags = [(2, tid["tag"]) for tid in rotor_dict.values()]
        stator_tags = [(2, tid["tag"]) for tid in stator_dict.values()]
        kwargs = dict(removeObject=True, removeTool=False)

        surf = gmsh_dict[idRotorAG]
        surf["tag"] = model.occ.cut([(2, surf["tag"])], rotor_tags, **kwargs)[0][0][1]

        surf = gmsh_dict[idStatorAG]
        surf["tag"] = model.occ.cut([(2, surf["tag"])], stator_tags, **kwargs)[0][0][1]

        factory.synchronize()

        # append air gap and sliding band to respective dict
        for idx, surf in gmsh_dict.items():
            if ROTOR_LAB_S in surf["label"]:
                rotor_dict.update({idx: gmsh_dict[idx]})
            if STATOR_LAB_S in surf["label"]:
                stator_dict.update({idx: gmsh_dict[idx]})

    ###################
    # Adding Airbox
    ###################
    gmsh_dict = {}

    if is_airbox and (not is_lam_only_R) and (not is_lam_only_S):
        ab_list = get_air_box(sym=sym, machine=machine)

        for surf in ab_list:
            nSurf += 1
            args = [boundary_prop, mesh_size_AB, user_mesh_dict]
            _draw_lines(model, surf, gmsh_dict, nSurf, *args)

        # draw surfaces (from lines) in gmsh and add surface tags to gmsh_dict
        for surf in gmsh_dict.values():
            _draw_surf(factory, surf)

    stator_dict.update(gmsh_dict)

    ######################
    ### Finalize Model ###
    ######################
    gmsh_dict = {}
    gmsh_dict.update(rotor_dict)
    gmsh_dict.update(stator_dict)

    # finally add physical groups
    for surf in gmsh_dict.values():
        model.addPhysicalGroup(2, [surf["tag"]], name=surf["label"])

    factory.synchronize()

    # while cutting surfaces, GMSH may have altered line tags
    # so we need to update our gmsh_dict
    _update_lines(output, model, gmsh_dict)

    # Set boundary conditions in gmsh lines
    boundary_list = list(set(boundary_prop.values()))
    for propname in boundary_list:
        bc_id = []
        for surf in gmsh_dict.values():
            for line in surf["lines"]:
                if line["bc_name"] == propname:
                    bc_id.extend([abs(line["tag"])])
        if bc_id:
            model.addPhysicalGroup(1, bc_id, name=propname)
            factory.synchronize()

        print(propname)
        print(bc_id)

    # Set all line labels as physical groups
    if is_set_labels:
        groups = {}
        for surf in gmsh_dict.values():
            for lvalues in surf.values():
                if (
                    type(lvalues) is not dict
                    or "label" not in lvalues
                    or not lvalues["label"]
                ):
                    continue

                if lvalues["label"] not in groups.keys():
                    groups[lvalues["label"]] = []
                groups[lvalues["label"]].append(abs(lvalues["tag"]))

        for label, tags in groups.items():
            factory.synchronize()
            model.addPhysicalGroup(1, tags, name=label)

    # save mesh or geo file depending on file extension
    filename, file_extension = splitext(path_save)

    if file_extension == ".geo":
        gmsh.write(filename + ".geo_unrolled")
        replace(filename + ".geo_unrolled", filename + file_extension)
    else:
        gmsh.model.mesh.generate(2)
        gmsh.write(path_save)

    if is_run:
        gmsh.fltk.run()  # Uncomment to launch Gmsh GUI
    gmsh.finalize()
    return gmsh_dict


def _draw_lines(model, surf, gmsh_dict, nSurf, bnd_prop, mesh_size, user_mesh):
    """
    Draw surface lines of surface with index nSurf in gmsh
    and store needed information in gmsh_dict[nSurf].
    Further store surfaces in tool_dict that have to been cut from main surface.
    """

    gmsh_dict.update(
        {nSurf: {"tag": None, "label": short_label(surf.label), "lines": []}}
    )

    # comp. number of elements on the lines & override by user values in case
    mesh_dict = comp_gmsh_mesh_dict(
        surface=surf, element_size=mesh_size, user_mesh_dict=user_mesh
    )

    # if surface is a SurfRing, first draw inner lines as a cutting tool
    if isinstance(surf, SurfRing):
        tool = surf.in_surf
        tool_dict = {0: {"tag": None, "label": None, "refSurfId": nSurf, "lines": []}}
        draw_surf_line(tool, mesh_dict, bnd_prop, model, tool_dict, 0, mesh_size)
        gmsh_dict[nSurf].update({"cut": tool_dict})
        surf = surf.out_surf

    # draw the surface lines
    draw_surf_line(surf, mesh_dict, bnd_prop, model, gmsh_dict, nSurf, mesh_size)


def _draw_surf(factory, surf):
    """
    Draw surfaces (from lines) in gmsh and add surface tags to gmsh_dict.
    """
    # build a lineloop of the surfaces lines
    lines = []
    for line in surf["lines"]:
        lines.extend([line["tag"]])
    cloop = factory.addCurveLoop(lines)

    # add surface
    surf["tag"] = factory.addPlaneSurface([cloop])
    factory.synchronize()

    # draw cutting tool and cut surface if needed
    if "cut" in surf:
        tool = surf["cut"][0]
        _draw_surf(factory, tool)

        kwargs = dict(removeObject=True, removeTool=True)
        cut = factory.cut([(2, surf["tag"])], [(2, tool["tag"])], **kwargs)
        surf["tag"] = cut[0][0][1]
        factory.synchronize()

    return None


def _update_lines(output, model, gmsh_dict, tolerance=1e-9):
    """Update line tags in gmsh_dict if possible."""
    for surf in gmsh_dict.values():

        # get the list of lines for the surface
        new_lines = _get_surf_line_list(model, surf)
        old_lines = surf["lines"]

        # check the original list of lines if there is a corresponding line
        # in the new line list of the GMSH model
        for old_line in old_lines:
            # readability
            p1_old = old_line["begin"]["coord"]
            p2_old = old_line["end"]["coord"]
            com_old = old_line["COM"]  # center of mass
            typ_old = old_line["typ"]

            update_count = 0
            for new_line in new_lines:
                # skip if the new line was already used for an update
                if "skip" in new_line:
                    continue

                # readability
                p1_new = new_line["coord"][0]
                p2_new = new_line["coord"][1]
                com_new = new_line["COM"]  # center of mass
                typ_new = new_line["typ"]

                # some conditions to identify matching lines
                condcc = _norm(com_old, com_new) <= tolerance
                cond11 = _norm(p1_old, p1_new) <= tolerance
                cond22 = _norm(p2_old, p2_new) <= tolerance
                cond21 = _norm(p2_old, p1_new) <= tolerance
                cond12 = _norm(p1_old, p2_new) <= tolerance
                condtyp = typ_old == typ_new

                if condcc and ((cond11 and cond22) or (cond12 and cond21)):
                    # lines seem to match perfectly so update tag
                    # ... but we don't care for sign
                    old_line["tag"] = new_line["tag"]
                    new_line["skip"] = True  # TODO is it valid to skip?
                    update_count += 1

                elif (cond11 or cond12 or cond21 or cond22) and condtyp:
                    # at least one point of line match, do further checks
                    if typ_new == "Line":
                        # ckeck if new line is on old line
                        if _is_line_inside(p1_new, p2_new, p1_old, p2_old):
                            print()

                    print()

            if update_count == 0 and old_line["bc_name"]:
                output.get_logger().warning(
                    "draw_gmsh warning: _update_line_tags() "
                    + "found no line to update tag - "
                    + f"old line tag: {old_line['tag']} - surface: {surf['label']}"
                )


def _norm(p1, p2):
    """Calculate the norm of 3 points in n dimensions."""
    sqr = [(x1 - x2) ** 2 for x1, x2 in zip(p1, p2)]
    return sum(sqr) ** 1 / 2


def _get_surf_line_list(model, surf):
    # get the list of lines for the surface
    lines = []
    dimTags = model.getBoundary([(2, surf["tag"])])
    for dimTag in dimTags:
        dimTag = [dimTag[0], abs(dimTag[1])]
        typ = model.getType(dimTag[0], dimTag[1])
        bnd = model.getBoundary([dimTag])
        coord = [model.getValue(0, dTag[1], []) for dTag in bnd]
        COM = model.occ.getCenterOfMass(*dimTag)
        lines_dict = dict(tag=dimTag[1], typ=typ, coord=coord, COM=COM)
        lines.append(lines_dict)

    return lines


def _is_on_line(p1, p2, q2):
    """Check if point q lies on line segment 'pr'"""
    if (
        p1[0] <= max(p2[0], q2[0])
        and p1[0] >= min(p2[0], q2[0])
        and p1[1] <= max(p2[1], q2[1])
        and p1[1] >= min(p2[1], q2[1])
    ):
        return True
    return False


def _is_collinear(p0, p1, p2, tolerance=1e-9):
    x1, y1 = p1[0] - p0[0], p1[1] - p0[1]
    x2, y2 = p2[0] - p0[0], p2[1] - p0[1]
    cond = abs(x1 * y2 - x2 * y1) < tolerance
    return cond


def _is_line_inside(p1, q1, p2, q2):
    """Return True if points p2 and q2 are completely inside 'p1q1'"""
    # Check if p2 and q2 are collinear with p1q1 and lie on segment p1q1
    if _is_collinear(p1, q1, p2) == 0 and _is_collinear(p1, q1, q2) == 0:
        if _is_on_line(p2, p1, q1) and _is_on_line(q2, p1, q1):
            return True

    return False
