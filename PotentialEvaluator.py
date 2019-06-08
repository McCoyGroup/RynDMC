import numpy as np
from CPotentialLib import CPotential

##################################################################################################################
##
##                                                  PotentialEvaluator
##
##

class PotentialEvaluator:

    def __init__(self, func, mode="single"):

        self.f = func
        if type(func).__name__ == "PyCapsule": #poor man's type check...
            self.f = CPotential(func, mode = mode)
        self.mode = mode

    def call(self, atoms, coords):
        if self.mode == "single":
            return np.array([ self.f(atoms, coord) for coord in coords ])
        else:
            return self.f(atoms, coords)

    def __call__(self, *args, **kwargs):
        return self.call(*args, **kwargs)