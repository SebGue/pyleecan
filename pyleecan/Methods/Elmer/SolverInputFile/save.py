# -*- coding: utf-8 -*-


def save(self, sif_file):
    """ Save the Elmer Solver Input File (sif) to a file."""
    with open(sif_file, "wt") as f:
        self.write(f)
