# -*- coding: utf-8 -*-
from os import replace
from os.path import splitext
from ....Functions.labels import (
    SHAFT_LAB,
    decode_label,
    HOLEM_LAB_S,
    HOLEV_LAB_S,
    ROTOR_LAB_S,
    ROTOR_LAB,
    MAG_LAB,
    STATOR_LAB_S,
    LAM_LAB_S,
)

from ....Methods.Simulation.StructElmer import MASTER, SLAVE

from numpy import array, dot, linalg


def process_mesh(
    self, gmsh_dict, lam_file, mag_file, is_get_magnet=False, is_hole_air=True
):
    """Preprocess the GMSH model, i.e. remove unused parts, rename boundaries, ..."""
    # TODO utilize 'translation' dict

    gmsh = gmsh_dict["gmsh"]
    model = gmsh_dict["model"]
    factory = gmsh_dict["factory"]

    # remove unused model parts
    sta_label = (STATOR_LAB_S + "-0_").lower()  # Stator parts
    _remove_entities(model, factory, labels=[sta_label])

    # get group names
    grps = model.getPhysicalGroups(-1)
    grp_names = [model.getPhysicalName(*grp) for grp in grps]

    # get lists of some surfaces by name
    magnet_list = []
    for grp, name in zip(grps, grp_names):
        label_dict = decode_label(name)
        if HOLEM_LAB_S in label_dict["surf_type"]:
            entities = model.getEntitiesForPhysicalGroup(*grp)
            if grp[0] == 2:
                magnet_list.extend(entities.tolist())

    lam_list = []
    for grp, name in zip(grps, grp_names):
        label = decode_label(name)
        if ROTOR_LAB_S in label["lam_type"] and LAM_LAB_S in label["surf_type"]:
            entities = model.getEntitiesForPhysicalGroup(*grp)
            if grp[0] == 2:
                lam_list.extend(entities.tolist())

    lam_lines = []
    for lam in lam_list:
        lam_lines.extend(model.getBoundary([(2, lam)], oriented=False))
    lam_lines = list(set([lam[1] for lam in lam_lines]))  # unique

    hole_lines = []
    for line in lam_lines:
        names = _get_names_physical(gmsh, dimtag=[1, line])
        if any([HOLEV_LAB_S in decode_label(name)["surf_type"] for name in names]):
            hole_lines.append(line)

    # setup dict to store physical groups, key: group name, value: list of tags
    groups_dict = {}

    # get lines of magnets for processing their physical groups
    for id, magnet in enumerate(magnet_list):
        lines = model.getBoundary([(2, magnet)])
        CoM = model.occ.getCenterOfMass(2, magnet)
        affected = _centrifugal_edges(model, lines, CoM)

        # store new group names in 'groups_dict' to set it later
        for ii, line in enumerate(affected):
            name = "_".join([ROTOR_LAB, MAG_LAB, f"{id}", f"{ii}"])
            if line[1] in lam_lines:  # only lines with direct contact for now
                if name not in groups_dict.keys():
                    groups_dict[name] = []
                groups_dict[name].append(line)

        # store new magnet body name
        key = "_".join(["Magnet", str(id), "Body"])
        groups_dict[key] = [(2, magnet)]

    # store new lamination body name
    for id, lam in enumerate(lam_list):
        s = "Lamination_" + str(id) + "_Body"
        if s not in groups_dict:
            groups_dict[s] = []
        groups_dict[s].append((2, lam))

    # store hole if not air
    if not is_hole_air:
        pass  # TODO

    # add boundaries to keep to groups_dict
    keeper_list = [
        "ROTOR" + MASTER,
        "ROTOR" + SLAVE,
        "Rotor_Tangential_Bridge",
        "Rotor_Radial_Bridge",
        "ROTOR_BORE_CURVE",
    ]

    for line in lam_lines:
        names = _get_names_physical(gmsh, dimtag=[1, line])
        for key in keeper_list:
            for name in names:
                if key in name:
                    if name not in groups_dict:
                        groups_dict[name] = []
                    groups_dict[name].append((1, line))

    # # update group names
    # grps = model.getPhysicalGroups(-1)
    # grp_names = [model.getPhysicalName(*grp) for grp in grps]

    # # delete unused surfaces
    # RL = ROTOR_LAB_S + "-0_"
    # del_list = [SHAFT_LAB, RL + HOLEV_LAB_S]
    # if not is_get_magnet:
    #     del_list.append(RL + HOLEM_LAB_S)

    # if not is_get_lam:
    #     rot_lam = (RL + LAM_LAB_S).lower()  # Rotor Lamination
    #     del_list.append(rot_lam)

    # for grp, name in zip(grps, grp_names):
    #     if any([n in name.lower() for n in del_list]):
    #         entities = model.getEntitiesForPhysicalGroup(*grp).tolist()
    #         for entity in entities:
    #             if grp[0] == 2:
    #                 factory.remove([(2, entity)], recursive=False)
    # factory.synchronize()

    # meshing
    surf_list = model.getEntities(dim=2)
    for surf in surf_list:
        model.mesh.setRecombine(2, surf[1])

    # gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
    model.mesh.generate(2)
    # model.mesh.recombine()
    model.mesh.refine()

    # remove all 'old' pyhsical groups and set new physical group names
    model.removePhysicalGroups(dimTags=[])
    for name in grp_names:
        model.removePhysicalName(name)

    mesh_names = dict()
    # magnets
    if is_get_magnet:
        mesh_names["Magnets"] = []
        items = [(key, item) for (key, item) in groups_dict.items() if "Magnet" in key]
        for key, values in items:
            for dim in (1, 2):
                tags = [abs(tag) for d, tag in values if d == dim]
                if not tags:
                    continue
                name = f"{key}_Master" if dim == 1 else key
                model.addPhysicalGroup(dim, tags, name=name)
                mesh_names["Magnets"].append(name)

        # save mesh or geo file depending on file extension
        filename, file_extension = splitext(mag_file)
        filename += "_magnets"

        if file_extension == ".geo":
            gmsh.write(filename + ".geo_unrolled")
            replace(filename + ".geo_unrolled", filename + file_extension)
        else:
            gmsh.write(mag_file)

    #
    model.removePhysicalGroups(dimTags=[])
    for name in grp_names:
        model.removePhysicalName(name)

    # lamination
    mesh_names["Lamination"] = []
    for key, values in groups_dict.items():
        for dim in [1, 2]:
            tags = [abs(dimtag[1]) for dimtag in values if dimtag[0] == dim]
            if not tags:
                continue
            if "Magnet" in key and is_get_magnet and dim == 1:
                key = key + "_Slave" if dim == 1 else key
                model.addPhysicalGroup(dim, tags, name=key)
                mesh_names["Lamination"].append(key)
            elif not "Magnet" in key:
                model.addPhysicalGroup(dim, tags, name=key)
                mesh_names["Lamination"].append(key)

    # save mesh or geo file depending on file extension
    filename, file_extension = splitext(lam_file)

    if file_extension == ".geo":
        gmsh.write(filename + ".geo_unrolled")
        replace(filename + ".geo_unrolled", filename + file_extension)
    else:
        gmsh.write(lam_file)

    return mesh_names


def _remove_entities(model, factory, labels):
    """
    Remove model entities that have one of the given labels in their physical
    groups names.

    Parameters
    ----------

    gmsh :
        gmsh object

    labels : list
        list of labels

    """
    # TODO add test that entities are not in surf or part of other 'keeper'entities
    # get all group names
    grps = model.getPhysicalGroups(-1)
    grp_names = [model.getPhysicalName(*grp) for grp in grps]

    # get entities that will be removed
    pt_list = []
    line_list = []
    surf_list = []
    for grp, name in zip(grps, grp_names):
        if any([label in name.lower() for label in labels]):
            dim = grp[0]
            entities = model.getEntitiesForPhysicalGroup(dim, grp[1])
            if dim == 0:
                pt_list.extend(entities.tolist())
            elif dim == 1:
                line_list.extend(entities.tolist())
            elif dim == 2:
                surf_list.extend(entities.tolist())

    # get lines of surfaces
    for surf in surf_list:
        lines = model.getBoundary([(2, surf)])  # TODO check if holes are included
        for line in lines:
            line_list.append(line[1])

    # get points of lines
    for line in line_list:
        pts = model.getBoundary([(1, line)])
        for pt in pts:
            pt_list.append(pt[1])

    # get unique list of entities to remove
    line_list = list(set(line_list))
    pt_list = list(set(pt_list))

    # delete unused entities
    for surf in surf_list:
        # model.removeEntities((2, surf), recursive=False)
        factory.remove([(2, surf)], recursive=False)

    for line in line_list:
        # model.removeEntities((1, line), recursive=False)
        factory.remove([(1, line)], recursive=False)

    for pt in pt_list:
        # model.removeEntities((0, pt), recursive=False)
        factory.remove([(0, pt)], recursive=False)

    # synchronize to apply changes to model
    factory.synchronize()


def _get_names_physical(gmsh, dimtag):
    model = gmsh.model
    grp_tags = model.getPhysicalGroupsForEntity(*dimtag)
    names = [model.getPhysicalName(1, tag) for tag in grp_tags]

    return names


def _centrifugal_edges(model, lines, com):
    """
    Determine lines of a surface on which a centrifugal force acts when the rotation
    is around the origin (0,0).

    lines: list of lines (dimTags)
    com: center of mass (x,y)

    Returns: list of affected lines
    """
    origin = array([0.0, 0.0])
    com = com[0:2]
    force = com - origin  # direction of force

    affected = []

    for line in lines:
        line = (line[0], abs(line[1]))
        p0 = model.getValue(*line, [0])[0:2]
        p1 = model.getValue(*line, [1])[0:2]
        mid = model.getValue(*line, [0.5])[0:2]
        edge_vec = p1 - p0

        # calculate normal
        n = array([-edge_vec[1], edge_vec[0]])
        n /= linalg.norm(n)

        # ensure normal points outward
        if dot(n, mid - com) < 0:
            n = -n

        # if the force points in the direction of the normal, it acts on the line
        if dot(force, n) > 0:
            affected.append(line)
    return affected
