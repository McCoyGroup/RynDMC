from McUtils.Data import UnitsData, AtomData
from Coordinerds.CoordinateSystems import CoordinateSet
import numpy as np


##################################################################################################################
##
##                                                  WalkerSet
##
##
class WalkerSet:
    def __init__(self,
                 atoms,
                 coords,
                 *miscargs,
                 initial_walkers=5000,
                 masses = None,
                 sigmas = None,
                 timestep = 1, # I don't like having this here but I can't really help it...
                 _initialize = True
                 ):

        if _initialize:
            if len(miscargs) > 0:
                raise ValueError("{}: only the `atoms` and `coords` arguments should be passed as non-keyword args".format(
                    type(self).__name__
                ))

            atoms = [ AtomData[atom, "Symbol"] for atom in atoms ] # take whatever atom spec and symbolize it
            if masses is None:
                masses = np.array([ UnitsData.convert(AtomData[atom, "Mass"], "AtomicUnitOfMass") for atom in atoms ])

            if sigmas is None:
                sigmas = np.sqrt(masses)

            coords = CoordinateSet(coords) # by default CartesianCoordinates3D

            if isinstance(initial_walkers, int):
                #out of laziness we'll just duplicate our original a bunch
                initial_walkers =  coords*initial_walkers
            else:
                initial_walkers = CoordinateSet(initial_walkers)

        # I could add some fancy property tricks to prevent this stuff from being mangled, but it's not worth it...
        self.atoms = atoms
        self.masses = masses
        self.timestep = timestep
        self.base_coords = coords
        self.sigmas = sigmas
        self.coords = initial_walkers

    def get_displacements(self, n=1, sigmas = None):
        """Computes n random Gaussian displacements from the sigmas and the timestep

        :param n:
        :type n:
        :param sigmas:
        :type sigmas:
        :return:
        :rtype:
        """
        if sigmas is None:
            sigmas = self.sigmas
        return np.array([
            np.random.normal(0.0, sig, (len(self.coords), 3)) for sig in sigmas
        ]).T

    def get_displaced_coords(self, n=1):
        """Computes n new displaced coordinates based on the original ones from the sigmas

        :param n:
        :type n:
        :return:
        :rtype:
        """
        accum_disp = np.cumsum(self.get_displacements(n))
        return self.coords + accum_disp # hoping the broadcasting makes this work...
