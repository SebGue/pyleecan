# -*- coding: utf-8 -*-
from .loggers import init_default_log
import os
import platform
import subprocess

PACKAGE_NAME = "pyleecan"
# User folder (to store machine/materials/config)
if platform.system() == "Windows":
    USER_DIR = os.path.join(os.environ["APPDATA"], PACKAGE_NAME)
    USER_DIR = USER_DIR.replace("\\", "/")
else:
    USER_DIR = os.environ["HOME"] + "/.local/share/" + PACKAGE_NAME

__version__ = "1.5.2"

try:
    __commit__ = (
        subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=os.path.dirname(os.path.abspath(__file__))
        )
        .strip()
        .decode("utf-8")
    )
except Exception:
    __commit__ = ""

init_default_log()
