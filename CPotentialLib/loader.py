import os, sys

lib_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lib_dir)
try:
    from .CPotentialLib import *
except ImportError:
    from .src.setup import load
    if load() is None:
        raise
    from .CPotentialLib import *
finally:
    sys.path.pop(0)