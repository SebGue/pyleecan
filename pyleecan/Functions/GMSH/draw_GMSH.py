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

from numpy import pi

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

        # finally add physical groups TODO maybe only after cutting
        for surf in gmsh_dict.values():
            model.addPhysicalGroup(2, [surf["tag"]], name=surf["label"])

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

        # finally add physical groups TODO maybe only after cutting
        for surf in gmsh_dict.values():
            model.addPhysicalGroup(2, [surf["tag"]], name=surf["label"])

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

        # finally add physical groups TODO maybe only after cutting
        for surf in gmsh_dict.values():
            model.addPhysicalGroup(2, [surf["tag"]], name=surf["label"])

        factory.synchronize()

        # append air gap and sliding band to stator and rotor
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

        # finally add physical groups TODO maybe only after cutting
        for surf in gmsh_dict.values():
            model.addPhysicalGroup(2, [surf["tag"]], name=surf["label"])

        factory.synchronize()

        stator_dict.update(gmsh_dict)

    print()
    # if is_sliding_band and (not is_lam_only_R) and (not is_lam_only_S):
    #     if sym == 1:
    #
    #     else:
    #
    #         # Look at the lines in the resulting surface, then update the dictionary
    #         # with MASTER/SLAVE BC when line angles match symmetry angles
    #         # MASTER is x-axis and SLAVE is 2Pi/sym
    #         rotor_ag_after = model.getEntitiesForPhysicalGroup(2, pg)
    #         rotor_ag_new_lines = model.getBoundary([(2, rotor_ag_after)])
    #         nline = 0
    #         for type_entity_l, rotor_ag_line in rotor_ag_new_lines:
    #             rotor_ag_new_points = model.getBoundary(
    #                 [(type_entity_l, abs(rotor_ag_line))]
    #             )
    #             btag = rotor_ag_new_points[0][1]
    #             etag = rotor_ag_new_points[1][1]
    #             bxy = model.getValue(0, btag, [])
    #             exy = model.getValue(0, etag, [])
    #             exy_angle = cmath.phase(complex(exy[0], exy[1]))
    #             bxy_angle = cmath.phase(complex(bxy[0], bxy[1]))
    #             if exy_angle == bxy_angle and abs(exy_angle) < 1e-6:
    #                 b_name = boundary_prop[AS_BR_LAB]
    #                 l_name = AS_BR_LAB
    #             elif (exy_angle == bxy_angle) and (
    #                 abs(abs(exy_angle) - 2.0 * pi / sym) < 1e-6
    #             ):
    #                 b_name = boundary_prop[AS_BL_LAB]
    #                 l_name = AS_BL_LAB
    #             else:
    #                 b_name = None
    #                 l_name = None
    #             gmsh_dict[rotor_ag_before[1]].update(
    #                 {
    #                     "tag": rotor_ag_after[0],
    #                     "label": lab_int + "_" + AIRGAP_LAB + BOT_LAB,
    #                     nline: {
    #                         "tag": abs(rotor_ag_line),
    #                         "label": l_name,
    #                         "n_elements": None,
    #                         "bc_name": b_name,
    #                         "begin": {
    #                             "tag": btag,
    #                             "coord": complex(bxy[0], bxy[1]),
    #                         },
    #                         "end": {"tag": etag, "coord": complex(exy[0], exy[1])},
    #                         "arc_angle": None,
    #                         "line_angle": None,
    #                     },
    #                 }
    #             )
    #             nline = nline + 1

    #         name = lab_ext + "_" + AIRGAP_LAB + TOP_LAB
    #         if len(cut2[0]) > 1:
    #             # Remove extra surfaces
    #             model.occ.remove([cut2[0][0]])
    #             factory.synchronize()
    #             pg = model.addPhysicalGroup(2, [cut2[0][1][1]], name=name)
    #         else:
    #             factory.synchronize()
    #             pg = model.addPhysicalGroup(2, [cut2[0][0][1]], name=name)

    #         # Look at the lines in the resulting surface, then update the dictionary
    #         # with MASTER/SLAVE BC when line angles match symmetry angles
    #         # MASTER is x-axis and SLAVE is 2Pi/sym
    #         stator_ag_after = model.getEntitiesForPhysicalGroup(2, pg)
    #         stator_ag_new_lines = model.getBoundary([(2, stator_ag_after)])
    #         nline = 0
    #         for type_entity_l, stator_ag_line in stator_ag_new_lines:
    #             stator_ag_new_points = model.getBoundary(
    #                 [(type_entity_l, abs(stator_ag_line))]
    #             )
    #             btag = stator_ag_new_points[0][1]
    #             etag = stator_ag_new_points[1][1]
    #             bxy = model.getValue(0, btag, [])
    #             exy = model.getValue(0, etag, [])
    #             exy_angle = cmath.phase(complex(exy[0], exy[1]))
    #             bxy_angle = cmath.phase(complex(bxy[0], bxy[1]))
    #             if exy_angle == bxy_angle and abs(exy_angle) < 1e-6:
    #                 b_name = boundary_prop[AS_TR_LAB]
    #                 l_name = AS_TR_LAB
    #             elif (exy_angle == bxy_angle) and (
    #                 abs(abs(exy_angle) - 2.0 * pi / sym) < 1e-6
    #             ):
    #                 b_name = boundary_prop[AS_TL_LAB]
    #                 l_name = AS_TL_LAB
    #             else:
    #                 b_name = None
    #                 l_name = None
    #             gmsh_dict[stator_ag_before[1]].update(
    #                 {
    #                     "tag": stator_ag_after[0],
    #                     "label": lab_ext + "_" + AIRGAP_LAB + TOP_LAB,
    #                     nline: {
    #                         "tag": abs(stator_ag_line),
    #                         "label": l_name,
    #                         "n_elements": None,
    #                         "bc_name": b_name,
    #                         "begin": {
    #                             "tag": btag,
    #                             "coord": complex(bxy[0], bxy[1]),
    #                         },
    #                         "end": {"tag": etag, "coord": complex(exy[0], exy[1])},
    #                         "arc_angle": None,
    #                         "line_angle": None,
    #                     },
    #                 }
    #             )
    #             nline = nline + 1

    # Set boundary conditions in gmsh lines
    boundary_list = list(set(boundary_prop.values()))
    for propname in boundary_list:
        bc_id = []
        for s_data in gmsh_dict.values():
            for lvalues in s_data.values():
                if isinstance(lvalues, dict) and lvalues["bc_name"] == propname:
                    bc_id.extend([abs(lvalues["tag"])])
        if bc_id:
            factory.synchronize()
            model.addPhysicalGroup(1, bc_id, name=propname)

    # Set all line labels as physical groups
    if is_set_labels:
        groups = {}
        for s_data in gmsh_dict.values():
            for lvalues in s_data.values():
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
