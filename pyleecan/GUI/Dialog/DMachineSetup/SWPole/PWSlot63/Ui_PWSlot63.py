# -*- coding: utf-8 -*-

# File generated according to PWSlot63.ui
# WARNING! All changes made in this file will be lost!
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from ......GUI.Tools.FloatEdit import FloatEdit

from pyleecan.GUI.Resources import pyleecan_rc


class Ui_PWSlot63(object):
    def setupUi(self, PWSlot63):
        if not PWSlot63.objectName():
            PWSlot63.setObjectName("PWSlot63")
        PWSlot63.resize(1084, 543)
        PWSlot63.setMinimumSize(QSize(630, 470))
        PWSlot63.setMaximumSize(QSize(16777215, 16777215))
        self.horizontalLayout = QHBoxLayout(PWSlot63)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.img_slot = QLabel(PWSlot63)
        self.img_slot.setObjectName("img_slot")
        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.img_slot.sizePolicy().hasHeightForWidth())
        self.img_slot.setSizePolicy(sizePolicy)
        self.img_slot.setMaximumSize(QSize(16777215, 16777215))
        self.img_slot.setPixmap(
            QPixmap(":/images/images/MachineSetup/WSlot/SlotW63_wind_int_rotor.png")
        )
        self.img_slot.setScaledContents(False)
        self.img_slot.setAlignment(Qt.AlignCenter)

        self.verticalLayout_2.addWidget(self.img_slot)

        self.txt_constraint = QTextEdit(PWSlot63)
        self.txt_constraint.setObjectName("txt_constraint")
        sizePolicy1 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(
            self.txt_constraint.sizePolicy().hasHeightForWidth()
        )
        self.txt_constraint.setSizePolicy(sizePolicy1)
        self.txt_constraint.setMaximumSize(QSize(16777215, 60))
        self.txt_constraint.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.txt_constraint.setTextInteractionFlags(
            Qt.TextSelectableByKeyboard | Qt.TextSelectableByMouse
        )

        self.verticalLayout_2.addWidget(self.txt_constraint)

        self.horizontalLayout.addLayout(self.verticalLayout_2)

        self.scrollArea = QScrollArea(PWSlot63)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollArea.setMinimumSize(QSize(270, 0))
        self.scrollArea.setMaximumSize(QSize(270, 16777215))
        self.scrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 268, 523))
        self.verticalLayout = QVBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.unit_W2 = QLabel(self.scrollAreaWidgetContents)
        self.unit_W2.setObjectName("unit_W2")

        self.gridLayout.addWidget(self.unit_W2, 2, 2, 1, 1)

        self.lf_W1 = FloatEdit(self.scrollAreaWidgetContents)
        self.lf_W1.setObjectName("lf_W1")

        self.gridLayout.addWidget(self.lf_W1, 1, 1, 1, 1)

        self.unit_W1 = QLabel(self.scrollAreaWidgetContents)
        self.unit_W1.setObjectName("unit_W1")

        self.gridLayout.addWidget(self.unit_W1, 1, 2, 1, 1)

        self.lf_H0 = FloatEdit(self.scrollAreaWidgetContents)
        self.lf_H0.setObjectName("lf_H0")

        self.gridLayout.addWidget(self.lf_H0, 3, 1, 1, 1)

        self.in_W2 = QLabel(self.scrollAreaWidgetContents)
        self.in_W2.setObjectName("in_W2")

        self.gridLayout.addWidget(self.in_W2, 2, 0, 1, 1)

        self.lf_H1 = FloatEdit(self.scrollAreaWidgetContents)
        self.lf_H1.setObjectName("lf_H1")

        self.gridLayout.addWidget(self.lf_H1, 4, 1, 1, 1)

        self.unit_H0 = QLabel(self.scrollAreaWidgetContents)
        self.unit_H0.setObjectName("unit_H0")

        self.gridLayout.addWidget(self.unit_H0, 3, 2, 1, 1)

        self.lf_W2 = FloatEdit(self.scrollAreaWidgetContents)
        self.lf_W2.setObjectName("lf_W2")

        self.gridLayout.addWidget(self.lf_W2, 2, 1, 1, 1)

        self.in_H2 = QLabel(self.scrollAreaWidgetContents)
        self.in_H2.setObjectName("in_H2")

        self.gridLayout.addWidget(self.in_H2, 5, 0, 1, 1)

        self.in_W1 = QLabel(self.scrollAreaWidgetContents)
        self.in_W1.setObjectName("in_W1")

        self.gridLayout.addWidget(self.in_W1, 1, 0, 1, 1)

        self.unit_W0 = QLabel(self.scrollAreaWidgetContents)
        self.unit_W0.setObjectName("unit_W0")

        self.gridLayout.addWidget(self.unit_W0, 0, 2, 1, 1)

        self.lf_H2 = FloatEdit(self.scrollAreaWidgetContents)
        self.lf_H2.setObjectName("lf_H2")

        self.gridLayout.addWidget(self.lf_H2, 5, 1, 1, 1)

        self.lf_W0 = FloatEdit(self.scrollAreaWidgetContents)
        self.lf_W0.setObjectName("lf_W0")

        self.gridLayout.addWidget(self.lf_W0, 0, 1, 1, 1)

        self.in_H0 = QLabel(self.scrollAreaWidgetContents)
        self.in_H0.setObjectName("in_H0")

        self.gridLayout.addWidget(self.in_H0, 3, 0, 1, 1)

        self.in_W0 = QLabel(self.scrollAreaWidgetContents)
        self.in_W0.setObjectName("in_W0")

        self.gridLayout.addWidget(self.in_W0, 0, 0, 1, 1)

        self.unit_H2 = QLabel(self.scrollAreaWidgetContents)
        self.unit_H2.setObjectName("unit_H2")

        self.gridLayout.addWidget(self.unit_H2, 5, 2, 1, 1)

        self.in_H1 = QLabel(self.scrollAreaWidgetContents)
        self.in_H1.setObjectName("in_H1")

        self.gridLayout.addWidget(self.in_H1, 4, 0, 1, 1)

        self.c_H1_unit = QComboBox(self.scrollAreaWidgetContents)
        self.c_H1_unit.addItem("")
        self.c_H1_unit.addItem("")
        self.c_H1_unit.setObjectName("c_H1_unit")

        self.gridLayout.addWidget(self.c_H1_unit, 4, 2, 1, 1)

        self.verticalLayout.addLayout(self.gridLayout)

        self.verticalSpacer = QSpacerItem(
            20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding
        )

        self.verticalLayout.addItem(self.verticalSpacer)

        self.g_output = QGroupBox(self.scrollAreaWidgetContents)
        self.g_output.setObjectName("g_output")
        self.g_output.setMinimumSize(QSize(200, 0))
        self.verticalLayout_3 = QVBoxLayout(self.g_output)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.out_Wlam = QLabel(self.g_output)
        self.out_Wlam.setObjectName("out_Wlam")

        self.verticalLayout_3.addWidget(self.out_Wlam)

        self.out_slot_height = QLabel(self.g_output)
        self.out_slot_height.setObjectName("out_slot_height")

        self.verticalLayout_3.addWidget(self.out_slot_height)

        self.out_yoke_height = QLabel(self.g_output)
        self.out_yoke_height.setObjectName("out_yoke_height")

        self.verticalLayout_3.addWidget(self.out_yoke_height)

        self.out_wind_surface = QLabel(self.g_output)
        self.out_wind_surface.setObjectName("out_wind_surface")

        self.verticalLayout_3.addWidget(self.out_wind_surface)

        self.out_tot_surface = QLabel(self.g_output)
        self.out_tot_surface.setObjectName("out_tot_surface")

        self.verticalLayout_3.addWidget(self.out_tot_surface)

        self.out_op_angle = QLabel(self.g_output)
        self.out_op_angle.setObjectName("out_op_angle")

        self.verticalLayout_3.addWidget(self.out_op_angle)

        self.out_tooth_width = QLabel(self.g_output)
        self.out_tooth_width.setObjectName("out_tooth_width")

        self.verticalLayout_3.addWidget(self.out_tooth_width)

        self.verticalLayout.addWidget(self.g_output)

        self.scrollArea.setWidget(self.scrollAreaWidgetContents)

        self.horizontalLayout.addWidget(self.scrollArea)

        QWidget.setTabOrder(self.lf_W0, self.lf_W1)
        QWidget.setTabOrder(self.lf_W1, self.lf_W2)
        QWidget.setTabOrder(self.lf_W2, self.lf_H0)
        QWidget.setTabOrder(self.lf_H0, self.lf_H1)
        QWidget.setTabOrder(self.lf_H1, self.lf_H2)
        QWidget.setTabOrder(self.lf_H2, self.txt_constraint)
        QWidget.setTabOrder(self.txt_constraint, self.scrollArea)

        self.retranslateUi(PWSlot63)

        QMetaObject.connectSlotsByName(PWSlot63)

    # setupUi

    def retranslateUi(self, PWSlot63):
        PWSlot63.setWindowTitle(QCoreApplication.translate("PWSlot63", "Form", None))
        self.img_slot.setText("")
        self.txt_constraint.setHtml(
            QCoreApplication.translate(
                "PWSlot63",
                '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN" "http://www.w3.org/TR/REC-html40/strict.dtd">\n'
                '<html><head><meta name="qrichtext" content="1" /><style type="text/css">\n'
                "p, li { white-space: pre-wrap; }\n"
                "</style></head><body style=\" font-family:'DejaVu Sans'; font-size:8.15094pt; font-weight:400; font-style:normal;\">\n"
                '<p align="center" style=" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;"><span style=" font-family:\'DejaVu Sans\'; font-size:10pt; font-weight:600; text-decoration: underline;">Constraints :</span></p>\n'
                '<p align="center" style=" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;"><span style=" font-family:\'DejaVu Sans\'; font-size:10pt;">W1 &gt; W0</span></p></body></html>',
                None,
            )
        )
        self.unit_W2.setText(QCoreApplication.translate("PWSlot63", "m", None))
        self.unit_W1.setText(QCoreApplication.translate("PWSlot63", "m", None))
        self.in_W2.setText(QCoreApplication.translate("PWSlot63", "W2 :", None))
        self.unit_H0.setText(QCoreApplication.translate("PWSlot63", "m", None))
        self.in_H2.setText(QCoreApplication.translate("PWSlot63", "H2 :", None))
        self.in_W1.setText(QCoreApplication.translate("PWSlot63", "W1 :", None))
        self.unit_W0.setText(QCoreApplication.translate("PWSlot63", "m", None))
        self.in_H0.setText(QCoreApplication.translate("PWSlot63", "H0 :", None))
        self.in_W0.setText(QCoreApplication.translate("PWSlot63", "W0 :", None))
        self.unit_H2.setText(QCoreApplication.translate("PWSlot63", "m", None))
        self.in_H1.setText(QCoreApplication.translate("PWSlot63", "H1 :", None))
        self.c_H1_unit.setItemText(
            0, QCoreApplication.translate("PWSlot63", "[rad]", None)
        )
        self.c_H1_unit.setItemText(
            1, QCoreApplication.translate("PWSlot63", "[\u00b0]", None)
        )

        self.g_output.setTitle(QCoreApplication.translate("PWSlot63", "Output", None))
        self.out_Wlam.setText(
            QCoreApplication.translate("PWSlot63", "Lamination width : ?", None)
        )
        self.out_slot_height.setText(
            QCoreApplication.translate("PWSlot63", "Slot height : ?", None)
        )
        self.out_yoke_height.setText(
            QCoreApplication.translate("PWSlot63", "Yoke height : ?", None)
        )
        self.out_wind_surface.setText(
            QCoreApplication.translate("PWSlot63", "Winding surface : ?", None)
        )
        self.out_tot_surface.setText(
            QCoreApplication.translate("PWSlot63", "Total surface : ?", None)
        )
        self.out_op_angle.setText(
            QCoreApplication.translate("PWSlot63", "Opening angle : ?", None)
        )
        self.out_tooth_width.setText(
            QCoreApplication.translate("PWSlot63", "Tooth average width : ?", None)
        )

    # retranslateUi
