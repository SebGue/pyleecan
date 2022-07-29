# Load the machine
from os.path import join
from pyleecan.Classes.MagNetwork import MagNetwork
from pyleecan.Classes.OPdq import OPdq
from pyleecan.Functions.load import load
from pyleecan.definitions import DATA_DIR
import matplotlib.pyplot as plt
import pytest
from os.path import join

from numpy import ones, pi, array, linspace, cos, sqrt, zeros

from pyleecan.Classes.Simu1 import Simu1
from pyleecan.Classes.InputCurrent import InputCurrent
from pyleecan.Classes.MagFEMM import MagFEMM
from pyleecan.Classes.MagFEMM import MagFEMM
from pyleecan.Classes.ForceMT import ForceMT
from pyleecan.Functions.Plot import dict_2D

is_show_fig = True
