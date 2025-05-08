# -*- coding: utf-8 -*-
"""File generated according to PHoleM64/gen_list.json
WARNING! All changes made in this file will be lost!
"""
from pyleecan.GUI.Dialog.DMachineSetup.SMHoleMag.PHoleM64.Ui_PHoleM64 import Ui_PHoleM64


class Gen_PHoleM64(Ui_PHoleM64):
    def setupUi(self, PHoleM64):
        """Abstract class to update the widget according to the csv doc"""
        Ui_PHoleM64.setupUi(self, PHoleM64)
        # Setup of in_R0
        txt = self.tr("""Corner radii""")
        self.in_R0.setWhatsThis(txt)
        self.in_R0.setToolTip(txt)

        # Setup of lf_R0
        self.lf_R0.validator().setBottom(0)
        txt = self.tr("""Corner radii""")
        self.lf_R0.setWhatsThis(txt)
        self.lf_R0.setToolTip(txt)

        # Setup of in_W1
        txt = self.tr("""Magnet width""")
        self.in_W1.setWhatsThis(txt)
        self.in_W1.setToolTip(txt)

        # Setup of lf_W1
        self.lf_W1.validator().setBottom(0)
        txt = self.tr("""Magnet width""")
        self.lf_W1.setWhatsThis(txt)
        self.lf_W1.setToolTip(txt)

        # Setup of in_W2
        txt = self.tr("""Hole width""")
        self.in_W2.setWhatsThis(txt)
        self.in_W2.setToolTip(txt)

        # Setup of lf_W2
        self.lf_W2.validator().setBottom(0)
        txt = self.tr("""Hole width""")
        self.lf_W2.setWhatsThis(txt)
        self.lf_W2.setToolTip(txt)

        # Setup of in_W3
        txt = self.tr("""Tooth width""")
        self.in_W3.setWhatsThis(txt)
        self.in_W3.setToolTip(txt)

        # Setup of lf_W3
        self.lf_W3.validator().setBottom(0)
        txt = self.tr("""Tooth width""")
        self.lf_W3.setWhatsThis(txt)
        self.lf_W3.setToolTip(txt)

        # Setup of in_H0
        txt = self.tr("""Distance from bore""")
        self.in_H0.setWhatsThis(txt)
        self.in_H0.setToolTip(txt)

        # Setup of lf_H0
        self.lf_H0.validator().setBottom(0)
        txt = self.tr("""Distance from bore""")
        self.lf_H0.setWhatsThis(txt)
        self.lf_H0.setToolTip(txt)

        # Setup of in_H2
        txt = self.tr("""Magnet height""")
        self.in_H2.setWhatsThis(txt)
        self.in_H2.setToolTip(txt)

        # Setup of lf_H2
        self.lf_H2.validator().setBottom(0)
        txt = self.tr("""Magnet height""")
        self.lf_H2.setWhatsThis(txt)
        self.lf_H2.setToolTip(txt)
