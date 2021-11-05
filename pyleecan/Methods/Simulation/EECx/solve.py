# -*- coding: utf-8 -*-


def solve(self, *args):
    """EECx abstract class solve method won't do anything."""
    self.get_logger().warning(
        "EECx.solve(): "
        + "This is the method of an abstract class and hence will do nothing."
    )
