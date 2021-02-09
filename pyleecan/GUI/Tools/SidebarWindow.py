import sys
from PySide2 import QtWidgets, QtGui, QtCore

"""
def _setColor(widget):
    widget.setAutoFillBackground(True)
    p = widget.palette()
    p.setColor(widget.backgroundRole(), QtCore.Qt.red)
    widget.setPalette(p)
"""
PriColor1 = (54, 116, 255)
PriColor2 = (68, 138, 255)
PriColor3 = (131, 185, 255)

SecColor1 = (79, 91, 98)
SecColor3 = (49, 54, 59)
SecColor2 = (35, 38, 41)

PriTextColor1 = (0, 0, 0)
PriTextColor2 = (20, 20, 20)

SecTextColor1 = (255, 255, 255)
SecTextColor2 = (180, 180, 180)



STYLE_nav_panel = f"""
    QFrame {{
        background-color: rgb{SecColor2};
        margin: 0px;
        padding: 0 0 0 0;
    }}
    QPushButton {{
        margin: 0px;
        border: none;
        color: rgb{SecTextColor2};
        font: bold;
        font-size: 20px;
        background-color: rgb{SecColor2};
        border-left: 5px solid rgb{SecColor3};
    }}
    QPushButton:checked {{
        border-left: 5px solid rgb{PriColor2};
        background-color: rgb{SecColor3};
    }}
    QPushButton:hover {{
        color: rgb{SecTextColor1};
    }}
    QPushButton:hover:pressed {{
        border-left: 5px solid rgb{PriColor2};
        background-color: rgb{PriColor2};
    }}
"""

STYLE_io_stack = f"""
    QPushButton {{
        padding: 0 10 0 10;
        color: rgb{SecTextColor1};
    	border: 2px solid rgb{PriColor1};
    	border-radius: 5px;	
    	background-color: rgb{PriColor1};
        min-width: 80px;
        min-height: 30px;
    }}
    QPushButton:hover {{
    	background-color: rgb{PriColor2};
    }}
    QPushButton:pressed {{
    	background-color: rgb{PriColor3};
    }}
    QPushButton:disabled {{
    	color: rgb{SecTextColor2};
    }}
    
    QLabel {{
    	color: rgb{SecTextColor1};
    }}

    QLineEdit {{
        padding: 0 10 0 10;
        color: rgb{SecTextColor1};
    	border: 2px solid rgb{SecColor1};
    	border-radius: 5px;	
        min-height: 30px;
    }}
    QLineEdit:focus {{
        border: 2px solid rgb{PriColor2};
    }}

    QTextEdit {{
        padding: 0 10 0 10;
        color: rgb{SecTextColor1};
    	border: 2px solid rgb{SecColor1};
    	border-radius: 5px;	
        min-height: 30px;
    }}
    QTextEdit:focus {{
        border: 2px solid rgb{PriColor2};
    }}

    QListWidget {{
        border: 2px solid rgb{SecColor1};
        color: rgb{SecTextColor1};
    }}
    QListWidget:focus {{
        border: 2px solid rgb{PriColor2};
    }}
    QListWidget::item {{
        color: rgb{SecTextColor1};
    }}
    QListWidget::item:hover {{
        background-color: rgb{SecColor1};
    }}
    QListWidget::item:selected {{
        background-color: rgb{PriColor1};
    }}
    QListWidget::item:disabled {{
        color: rgb{SecTextColor2};
    }}


    QComboBox {{
        padding: 0 10 0 10;
        color: rgb{SecTextColor1};
    	border: 2px solid rgb{PriColor1};
    	border-radius: 5px;	
    	background-color: rgb{PriColor1};
        min-height: 30px;
    }}
    QComboBox:disabled {{
    	color: rgb{SecTextColor2};
    }}
    QComboBox:hover {{
    	background-color: rgb{PriColor2};
    }}
    QComboBox QAbstractItemView {{
        border: 2px solid rgb{PriColor1};
        border-radius: 5px;
        selection-background-color: rgb{PriColor1};
        color: rgb{SecTextColor1};
    }}

    QSpinBox {{
        padding: 0 10 0 10;
        color: rgb{SecTextColor1};
    	border: 2px solid rgb{SecColor1};
    	border-radius: 5px;	
        min-height: 30px;
    }}
    QSpinBox:focus {{
        border: 2px solid rgb{PriColor2};
    }}

    QCheckBox {{
        color: rgb{SecTextColor1};
    }}
    QCheckBox::indicator {{
        border: 2px solid rgb{SecColor1};
        width: 20px;
        height: 20px;
        border-radius: 5px;
    }}
    QCheckBox::indicator:hover {{
        border: 2px solid rgb{PriColor1};
    }}
    QCheckBox::indicator:checked {{
        background: 2px solid rgb{PriColor1};
    }}
    

    QGroupBox {{
        border: 2px solid rgb{SecColor1};
        border-radius: 5px;
        margin-top: 10px;
        padding: 2px 0px 0px 0px;
    }}
    QGroupBox::title {{
        color: rgb{SecTextColor1};
        subcontrol-origin: margin;
        subcontrol-position: top middle;
        padding: -10px 0px 0px 0px;
        }}
    QGroupBox::indicator {{
        border: 2px solid rgb{SecColor1};
        width: 15px;
        height: 15px;
        border-radius: 5px;
        background-color: rgb{SecColor3};
    }}
    QGroupBox::indicator:hover {{
        background-color: rgb{PriColor2};
    }}
    QGroupBox::indicator:checked {{
        background-color: rgb{PriColor1};
    }}
"""


class SidebarWindow(QtWidgets.QMainWindow):
    def __init__(self):
        # === App-Init ===
        super(SidebarWindow, self).__init__()
        self._title = "Pyleecan"
        self.setWindowTitle(self._title)
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)

        # === Main Widgets ===
        # Navigation Panel with Button Group
        self.nav_panel = QtWidgets.QFrame()
        self.nav_panel.setStyleSheet(STYLE_nav_panel)

        self.nav_btn_grp = QtWidgets.QButtonGroup()
        self.nav_btn_grp.setExclusive(True)
        self.nav_btn_grp.buttonClicked[int].connect(self.switch_stack)
        self.btn_grp_fct = []

        self.nav_layout = QtWidgets.QVBoxLayout(self.nav_panel)
        self.nav_layout.setSpacing(0)
        self.nav_layout.setContentsMargins(0, 0, 0, 0)
        self.nav_layout.addStretch(1)  # add stretch first

        # Sub Window Stack
        self.io_stack = QtWidgets.QStackedWidget(self)
        self.io_stack.setStyleSheet(STYLE_io_stack)

        # === Main Layout ===
        main_layout = QtWidgets.QHBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(self.nav_panel)
        main_layout.addWidget(self.io_stack)

        self._main.setLayout(main_layout)
        self._main.setStyleSheet(f"background-color: rgb{SecColor3}; font-size: 18px;")

        self.show()
        self.centerOnScreen()

    def close_application(self):
        sys.exit()

    def switch_stack(self, btn):
        # print('Button Nbr. %2d pressed' % btn)
        self.io_stack.setCurrentIndex(btn)  # set stack
        if self.btn_grp_fct[btn] is not None:  # execute user function
            self.btn_grp_fct[btn]()

    def centerOnScreen(self):
        """centerOnScreen() - Centers the window on the screen."""
        resolution = QtWidgets.QDesktopWidget().screenGeometry()
        frame = self.frameSize()
        self.move(
            (resolution.width() / 2) - (frame.width() / 2),
            (resolution.height() / 2) - (frame.height() / 2),
        )

    def addSubWindow(self, name, widget, btn_fct=None):
        """ add a new sub window to the stack including the coresponding button"""
        # Button
        btn = QtWidgets.QPushButton(name)
        btn.setFixedSize(120, 60)
        btn.setCheckable(True)

        self.nav_btn_grp.addButton(btn, self.io_stack.count())
        self.nav_layout.insertWidget(self.io_stack.count(), btn)
        self.btn_grp_fct.insert(self.io_stack.count(), btn_fct)

        # Stack
        self.io_stack.addWidget(widget)

        # enable first menu
        if self.io_stack.count() == 1:
            btn.setChecked(True)

    def eventFilter(self, obj, event):
        """
        Event Filter to disable 'Esc'-Key in a Widgets.
        To install eventFilter on a Widget:
            widget.installEventFilter(instance_of_main_window)
        """
        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() in (QtCore.Qt.Key_Escape,):
                return True
        return super(SidebarWindow, self).eventFilter(obj, event)

    def closeEvent(self, event):
        """Overload the methode to call DesignWidget.closeEvent """
        self.DesignWidget.closeEvent(event)
