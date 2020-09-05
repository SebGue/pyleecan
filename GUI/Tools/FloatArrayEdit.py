# -*- coding: utf-8 -*-
"""
@date Created on 2020-01-10 16:13:41
@author sebastian_g
@todo unittest it
"""
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import Qt

from pyleecan.GUI.Tools.FloatEdit import FloatEdit
from pyleecan.GUI import gui_option


class FloatArrayEdit(QtWidgets.QWidget):
    def __init__(self, unit="", *args, **kwargs):
        super(FloatArrayEdit, self).__init__(*args, **kwargs)

        # Set CSS
        self.setStyleSheet("QLineEdit { border: 0px }")

        # Create table
        self.tableWidget = QtWidgets.QTableWidget()
        self.tableWidget.setEditTriggers(QtWidgets.QTableWidget.NoEditTriggers)
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setRowCount(1)

        # Set the table headers
        self.tableWidget.setHorizontalHeaderLabels(["Header 1", "Header 2"])

        # Set the tooltips to headings
        self.tableWidget.horizontalHeaderItem(0).setToolTip("Column 1 ")
        self.tableWidget.horizontalHeaderItem(1).setToolTip("Column 2 ")

        # Set the alignment to the headers
        self.tableWidget.horizontalHeaderItem(0).setTextAlignment(Qt.AlignLeft)
        self.tableWidget.horizontalHeaderItem(1).setTextAlignment(Qt.AlignHCenter)

        # Set Data
        self.tableWidget.setCellWidget(0, 0, FloatEdit())
        self.tableWidget.setCellWidget(0, 1, FloatEdit())
        # self.tableWidget.move(0, 0)

        self.tableWidget.setFixedWidth(300)
        self.tableWidget.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.Stretch
        )

        # table selection change
        # self.tableWidget.doubleClicked.connect(self.on_click)

        # Button
        self.btn = QtWidgets.QPushButton("Print")
        self.btn.clicked.connect(self.on_click)

        # Layout
        self.setLayout(QtWidgets.QGridLayout())
        self.layout().addWidget(self.tableWidget, 0, 0)
        self.layout().addWidget(self.btn, 0, 1)
        # self.layout().setColumnStretch(0, 1)

    def on_click(self):
        print("Button Clicked\n")
        for currentQTableWidgetItem in self.tableWidget.selectedIndexes():
            row = currentQTableWidgetItem.row()
            col = currentQTableWidgetItem.column()
            val = self.tableWidget.cellWidget(row, col).value()
            print(row, col, val)
