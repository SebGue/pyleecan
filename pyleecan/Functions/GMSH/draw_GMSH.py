from ...Functions.labels import (
    short_label,
    LAM_LAB_S,
    ROTOR_LAB_S,
    STATOR_LAB_S,
    AIRGAP_LAB,
    LAM_LAB,
    STATOR_LAB,
    YOKE_LAB,
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

from numpy import pi


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
    is_mesh=True,
    is_finalize=True,
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
    is_mesh : bool
        True to run mesh generation and save mesh file except geo file is requested
    is_finalize : bool
        True to finalize model creation, otherwise model object is returned in output dict
        
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

    if user_mesh_dict is None:
        user_mesh_dict = {}

    # get machine
    machine = output.simu.machine

    # Default stator mesh element size
    mesh_size = 2 * pi * machine.rotor.Rext / 180
    mesh_size_min = mesh_size / 4
    mesh_size_max = mesh_size
    mesh_size_SB = 2.0 * pi * machine.rotor.Rext / 360.0  # Sliding Band
    mesh_size_AB = mesh_size / 2  # AirBox

    # For readibility
    model = gmsh.model
    factory = model.occ

    # Start a new model
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", int(False))
    gmsh.option.setNumber("Geometry.CopyMeshingMethod", 1)
    gmsh.option.setNumber("Geometry.PointNumbers", 0)
    gmsh.option.setNumber("Geometry.LineNumbers", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)
    gmsh.option.setNumber("Mesh.MeshSizeFactor", 3)
    model.add("Pyleecan")

    alpha = 0

    # remove VP0 BC if air box is used
    if is_airbox and (not is_lam_only_R) and (not is_lam_only_S):
        lab_ext = machine.get_lam_list()[1].get_label()
        boundary_prop = boundary_prop.copy()
        boundary_prop.pop(lab_ext + "_" + LAM_LAB + YOKE_LAB)

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
            args = [boundary_prop, mesh_size, user_mesh_dict]
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
            args = [boundary_prop, mesh_size, user_mesh_dict]
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
    rotor_air_dict = {}
    stator_air_dict = {}

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

        # store air gap and sliding band to respective dict
        for idx, surf in gmsh_dict.items():
            if ROTOR_LAB_S in surf["label"]:
                rotor_air_dict.update({idx: gmsh_dict[idx]})
            if STATOR_LAB_S in surf["label"]:
                stator_air_dict.update({idx: gmsh_dict[idx]})

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

    airbox_dict = gmsh_dict.copy()

    ######################
    ### Finalize Model ###
    ######################

    # dict of all surfaces
    gmsh_dict = {
        **rotor_dict,
        **rotor_air_dict,
        **stator_dict,
        **stator_air_dict,
        **airbox_dict,
    }

    # get all surfaces
    dimTagsRotor = []
    for surf in rotor_dict.values():
        dimTagsRotor.append((2, surf["tag"]))

    dimTagsStator = []
    for surf in stator_dict.values():
        dimTagsStator.append((2, surf["tag"]))

    dimTagsRotorAir = []
    for surf in rotor_air_dict.values():
        dimTagsRotorAir.append((2, surf["tag"]))

    dimTagsStatorAir = []
    for surf in stator_air_dict.values():
        dimTagsStatorAir.append((2, surf["tag"]))

    dimTagsAirBox = []
    for surf in airbox_dict.values():
        dimTagsAirBox.append((2, surf["tag"]))

    # make all interfaces coherent (except sliding band)
    factory.fragment(dimTagsRotor, dimTagsRotorAir)
    factory.fragment(dimTagsStator, dimTagsStatorAir)
    factory.fragment(dimTagsStator, dimTagsAirBox)
    factory.fragment(dimTagsRotor, dimTagsAirBox)
    factory.synchronize()

    # finally add surface physical groups
    for surf in gmsh_dict.values():
        model.addPhysicalGroup(2, [surf["tag"]], name=surf["label"])

    factory.synchronize()

    # set default boundaries on rotor and stator seperately to avoid issues
    rotor_combined_dict = {**rotor_dict, **rotor_air_dict}
    stator_combined_dict = {**stator_dict, **stator_air_dict}
    if machine.stator.is_outwards():
        stator_combined_dict.update(airbox_dict)
    else:
        rotor_combined_dict.update(airbox_dict)

    _set_default_boundaries(output, model, factory, rotor_combined_dict, boundary_prop)
    _set_default_boundaries(output, model, factory, stator_combined_dict, boundary_prop)

    # try to recover element sizes from original lines and set line element sizes
    _set_element_size(output, model, factory, gmsh_dict)

    # set all line labels as physical groups
    # TODO add 'not in boundary_list' since they have been set already
    if is_set_labels:
        pass  # TODO
        # groups = {}
        # for surf in gmsh_dict.values():
        #     for lvalues in surf.values():
        #         if (
        #             type(lvalues) is not dict
        #             or "label" not in lvalues
        #             or not lvalues["label"]
        #         ):
        #             continue

        #         if lvalues["label"] not in groups.keys():
        #             groups[lvalues["label"]] = []
        #         groups[lvalues["label"]].append(abs(lvalues["tag"]))

        # for label, tags in groups.items():
        #     factory.synchronize()
        #     model.addPhysicalGroup(1, tags, name=label)

    # save mesh or geo file depending on file extension
    filename, file_extension = splitext(path_save)

    if file_extension == ".geo":
        gmsh.write(filename + ".geo_unrolled")
        replace(filename + ".geo_unrolled", filename + file_extension)
    elif is_mesh:
        gmsh.model.mesh.generate(2)
        gmsh.write(path_save)

    if is_run:
        gmsh.fltk.run()

    if is_finalize:
        gmsh.finalize()
    else:
        gmsh_dict['gmsh'] = gmsh
        gmsh_dict['model'] = model
        gmsh_dict['factory'] = factory

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

    # common args for comp_gmsh_mesh_dict
    kwargs = dict(element_size=mesh_size, user_mesh_dict=user_mesh)

    # if surface is a SurfRing, first draw inner lines as a cutting tool
    if isinstance(surf, SurfRing):
        tool = surf.in_surf
        tool_dict = {0: {"tag": None, "label": None, "refSurfId": nSurf, "lines": []}}
        # comp. number of elements on the lines & override by user values in case
        mesh_dict = comp_gmsh_mesh_dict(surface=tool, **kwargs)
        # draw the surface lines
        draw_surf_line(tool, mesh_dict, bnd_prop, model, tool_dict, 0, mesh_size)
        gmsh_dict[nSurf].update({"cut": tool_dict})
        surf = surf.out_surf

    # comp. number of elements on the lines & override by user values in case
    mesh_dict = comp_gmsh_mesh_dict(surface=surf, **kwargs)

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


def _set_default_boundaries(output, model, factory, gmsh_dict, boundary_prop):
    # get all (orginal) lines and all new lines (without duplicates)
    # still line tags are not up to date due to cutting surfaces
    line_list = []
    for surf in gmsh_dict.values():
        line_list.extend(surf["lines"])
        if "cut" in surf:
            for tool in surf["cut"].values():
                line_list.extend(tool["lines"])

    new_line_list = []
    used_tags = set()
    for surf in gmsh_dict.values():
        for line in _get_surf_line_list(model, surf):
            tag = line["tag"]
            if tag not in used_tags:
                used_tags.add(tag)
                new_line_list.append(line)

    # set default boundary conditions in gmsh lines
    boundary_list = list(set(boundary_prop.values()))
    for propname in boundary_list:
        # get relevant lines
        bc_line_list = []
        for line in line_list:
            if line["bc_name"] == propname:
                bc_line_list.append(line.copy())

        # update line tag
        new_bc_lines = []
        for line in bc_line_list:
            lines = _get_updated_lines(output, line, new_line_list)
            new_bc_lines.extend(lines)

        if new_bc_lines:
            tags = list(set([line["tag"] for line in new_bc_lines]))  # remove dup.
            model.addPhysicalGroup(1, tags, name=propname)
            factory.synchronize()


def _get_updated_lines(output, line, new_lines, tolerance=1e-9):
    """Get the tags of actual lines in gmsh that are on 'line' or match 'line' completely."""

    # readability
    p1 = line["begin"]["coord"]
    p2 = line["end"]["coord"]
    com = line["COM"]  # center of mass
    typ = line["typ"]

    # loop all new lines to get original bc
    updated = []
    for new_line in new_lines:
        # readability
        p1_new = new_line["coord"][0]
        p2_new = new_line["coord"][1]
        com_new = new_line["COM"]  # center of mass
        typ_new = new_line["typ"]

        # some conditions to identify matching lines
        condcc = _norm(com, com_new) <= tolerance
        cond11 = _norm(p1, p1_new) <= tolerance
        cond22 = _norm(p2, p2_new) <= tolerance
        cond21 = _norm(p2, p1_new) <= tolerance
        cond12 = _norm(p1, p2_new) <= tolerance
        condtyp = typ == typ_new

        if condcc and ((cond11 and cond22) or (cond12 and cond21)):
            # lines seem to match perfectly so we don't need to search further
            return [new_line]

        if condtyp:
            # at least one point of line match, do further checks
            if typ_new == "Line":
                # ckeck if new line is on old line
                cond1 = _is_collinear(p1, p2, p1_new)
                cond2 = _is_collinear(p1, p2, p2_new)
                if cond1 and cond2:
                    cond3 = _is_on_line(p1_new, p1, p2)
                    cond4 = _is_on_line(p2_new, p1, p2)
                    if cond3 and cond4:
                        updated.append(new_line)
            elif typ_new == "Circle":
                pass  # TODO

    if not updated:
        output.get_logger().warning(
            "draw_gmsh() - Warning: "
            + "Found no corresponding line"
            + f" for {line['typ']} {line['tag']} with BC '{line['bc_name']}''."
        )

    return updated


def _set_element_size(output, model, factory, gmsh_dict):
    """Try to inherit element sizes for each surface line."""
    # get all (orginal) lines and all new lines (without duplicates)
    # still line tags are not up to date due to cutting surfaces
    size_dict = {}
    for surf in gmsh_dict.values():
        # set size seperately on each surface to avoid conflicts
        line_list = surf["lines"]
        if "cut" in surf:
            for tool in surf["cut"].values():
                line_list.extend(tool["lines"])

        # get the actual lines of the surface
        new_line_list = _get_surf_line_list(model, surf)

        # set element size to new line if possible
        for line in line_list:
            # skip if number of elements is not set
            if line["n_elements"] == 0:
                continue

            lines = _get_updated_lines(output, line, new_line_list)

            n_elem = line["n_elements"]
            if line["typ"] == "Line":
                len_old = _norm(line["begin"]["coord"], line["end"]["coord"])
                for line_new in lines:
                    tag = line_new["tag"]
                    len_new = _norm(line_new["coord"][0], line_new["coord"][1])
                    n_elem_new = int(max(round(n_elem / len_old * len_new), 1))
                    if tag not in size_dict or size_dict[tag] < n_elem_new:
                        model.mesh.setTransfiniteCurve(
                            tag, n_elem_new + 1, "Progression"
                        )
                        size_dict[tag] = n_elem_new
                    factory.synchronize()
                    # print(
                    #     f"Surface '{surf['label']}' Line {tag} set TransfiniteCurve {n_elem_new} "
                    # )
            else:
                for line_new in lines:
                    tag = line_new["tag"]
                    if tag not in size_dict or size_dict[tag] < n_elem:
                        model.mesh.setTransfiniteCurve(tag, n_elem + 1, "Progression")
                        size_dict[tag] = n_elem
                    factory.synchronize()


def _norm(p1, p2):
    """Calculate the norm of 3 points in n dimensions."""
    sqr = [(x1 - x2) ** 2 for x1, x2 in zip(p1, p2)]
    return sum(sqr) ** (1 / 2)


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


def _is_on_line(p1, p2, q2, tolerance=1e-9):
    """Check if point p1 lies on line segment 'p2q2'"""
    if (
        p1[0] <= (max(p2[0], q2[0]) + tolerance)
        and p1[0] >= (min(p2[0], q2[0]) - tolerance)
        and p1[1] <= (max(p2[1], q2[1]) + tolerance)
        and p1[1] >= (min(p2[1], q2[1]) - tolerance)
    ):
        return True
    return False


def _is_collinear(p0, p1, p2, tolerance=1e-9):
    x1, y1 = p1[0] - p0[0], p1[1] - p0[1]
    x2, y2 = p2[0] - p0[0], p2[1] - p0[1]
    cond = abs(x1 * y2 - x2 * y1) < tolerance
    return cond
