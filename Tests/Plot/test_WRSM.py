# -*- coding: utf-8 -*-
from os.path import join

from pyleecan.definitions import DATA_DIR
from pyleecan.Functions.load import load
from Tests import save_plot_path as save_path


def test_plot_WRSM():
    """Load and Plot the WRSM machine"""
    WRSM = load(join(DATA_DIR, "Machine", "WRSM_001.json"))
    WRSM.plot(is_show_fig=False, save_path=join(save_path, "test_WRSM.png"))


if __name__ == "__main__":
    test_plot_WRSM()
    print()
