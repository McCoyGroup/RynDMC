

class CPotential:
    def __init__(self, pointer, mode="single"):
        self.pointer = pointer # this should be a PyCapsule object
        self.mode = mode
    def call(self, atoms, coords):
        if self.mode == "single":
            from .loader import callPot
            return callPot(self.pointer, atoms, coords)
        elif self.mode == "mpi":
            from .loader import callMPIPot
            return callMPIPot(self.pointer, atoms, coords)
        else:
            from .loader import callPotVec
            return callPotVec(self.pointer, atoms, coords)
    def __call__(self, *args, **kwargs):
        return self.call(*args, **kwargs)