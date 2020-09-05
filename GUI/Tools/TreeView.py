# -*- coding: utf-8 -*-

from PyQt5 import QtCore, QtGui, QtWidgets
from os.path import join, isfile
from json import load as jload

from pyleecan.definitions import ROOT_DIR

def _get_prop(class_dict, cls_name, prop_type):
    """function to get a class property from a class_dict by name
    """
    if cls_name in class_dict.keys():
        if "properties" in class_dict[cls_name].keys():
            props = class_dict[cls_name]["properties"]
            for prop in props:
                if prop['name'] == prop_type:
                    return prop
        if "mother" in class_dict[cls_name].keys():
            if class_dict[cls_name]['mother'] != "":
                return _get_prop(class_dict, class_dict[cls_name]['mother'], prop_type)

    return None
    

# === TreeView ========================================================================
class TreeView(QtWidgets.QTreeView):
    def __init__(self):
        super(QtWidgets.QTreeView, self).__init__()

        # self.__typeWidgets = {}
        # self.__typeWrappers = {}

        model = QtGui.QStandardItemModel(0, 1)
        model.setHorizontalHeaderLabels(
            ["Name", "Value", "Unit", "Value Class", "Description"]
        )

        self.rootNode = model.invisibleRootItem()

        self.setModel(model)
        self.setColumnWidth(0, 150)

        self.setAlternatingRowColors(True)

        file_path = join(ROOT_DIR, "pyleecan", "Classes", "class_dict.json")
        # check if the file exist
        if not isfile(file_path):
            raise TreeViewMissingFileError(str(file_path) + " doesn't exist")

        # get the class dictionary
        with open(file_path, "r") as load_file:
            self.class_dict = jload(load_file)

    def generate(self, data, parent=None):
        if parent is None:
            parent = self.rootNode
            self.rootNode.removeRows(0, self.rootNode.rowCount())
        if hasattr(data, "as_dict"):
            self.generate_tree(data, parent)
        else:
            pass

    def generate_tree(self, data, parent):
        for attr in data.as_dict():
            if attr[0] != "_":
                _attr = getattr(data, attr)
                # --- pyleecan object attributes ---
                if hasattr(_attr, "as_dict"):
                    self.gen_branch(data, attr, parent)

                # --- float, int, str, none type attributes ---
                elif isinstance(_attr, (float, int, str, type(None))):
                    attribute = QtGui.QStandardItem(attr)
                    attribute.setEditable(False)
                    value = QtGui.QStandardItem(str(getattr(data, attr)))
                    value.setEditable(False)
                    # try to gather unit from class_dict, set unit "na" as default
                    cls_name = type(data).__name__
                    prop = _get_prop(self.class_dict, cls_name, attr)
                    unit_ = prop['unit'] if prop else 'na'
                    unit = QtGui.QStandardItem(unit_)

                    prop_type = QtGui.QStandardItem(type(_attr).__name__)
                    prop_doc = getattr(
                        type(data), attr
                    ).__doc__  # tc.__class__.prop.__doc__
                    prop_doc = QtGui.QStandardItem(prop_doc)
                    parent.appendRow([attribute, value, unit, prop_type, prop_doc])

                elif isinstance(_attr, list):
                    for elem in _attr:
                        pass

                # --- dict, ndarray type attributes ---
                else:
                    attribute = QtGui.QStandardItem(attr)
                    attribute.setEditable(False)
                    value = QtGui.QStandardItem('')
                    value.setEditable(False)
                    unit = QtGui.QStandardItem('')
                    prop_type = QtGui.QStandardItem(type(_attr).__name__)
                    prop_doc = QtGui.QStandardItem('not implemented yet')
                    parent.appendRow([attribute, value, unit, prop_type, prop_doc])

    def gen_branch(self, data, attr_name, parent):
        _attr = getattr(data, attr_name)
        branch = QtGui.QStandardItem(attr_name)
        branch.setEditable(False)
        value = QtGui.QStandardItem("")
        value.setEditable(False)
        # try to gather unit from class_dict, set unit "na" as default
        cls_name = type(data).__name__
        prop = _get_prop(self.class_dict, cls_name, attr_name)
        unit_ = prop['unit'] if prop else 'na'
        unit = QtGui.QStandardItem(unit_)

        prop_type = QtGui.QStandardItem(type(_attr).__name__)
        # class description rather than attribure description
        prop_doc = (
            getattr(type(data), attr_name).__doc__
            if _attr.__doc__ is None
            else _attr.__doc__
        )
        prop_doc = QtGui.QStandardItem(prop_doc)
        parent.appendRow([branch, value, unit, prop_type, prop_doc])
        self.generate(_attr, parent=branch)

class TreeViewMissingFileError(Exception):
    """ """

    pass