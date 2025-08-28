# -*- coding: utf-8 -*-

import os


def get_path_save(self, output):
    """Return the path to save the Elmer simulation

    Parameters
    ----------
    self : MagElmer
        a MagElmer object
    output : Output
        an Output object

        Returns
    -------
    save_path: str
        Path to save the Elmer simulation

    """
    path_res = output.get_path_result()

    save_dir = path_res
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    return save_dir
