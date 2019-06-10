from abc import *
import numpy as np
from .WalkerSet import WalkerSet
from .PotentialEvaluator import PotentialEvaluator


##################################################################################################################
##
##                                                  AbstractDMC
##
##
class AbstractDMC(metaclass=ABCMeta):
    """A completely abstract implementation of DMC"""

    def __init__(self, walker_set, potential, stack_len = None, log_file = None):
        if not isinstance(walker_set, WalkerSet):
            walker_set=WalkerSet(*walker_set)
        self.walkers = walker_set
        if not isinstance(potential, PotentialEvaluator):
            potential = PotentialEvaluator(potential)
        self.potential = potential
        if stack_len is None:
            self.stack = None
        else:
            from collections import deque
            self.stack = deque(maxlen=stack_len)
                # basically a circular buffer that we can feed stuff into in O(1) time
                # this out to get used for implementing descendant weighting in propogate
                # as memory can be cheaply tracked by feeding it into here

        self.log = log_file # we'll dump to binary here
        self.energies = [] # this'll fill out


    def diffuse(self, n=1):
        """Handles diffusion of walkers in DMC

        :param n: number of steps to generate diffusions for
        :type n:
        :return:
        :rtype:
        """
        newsies = self.walkers.get_displaced_coords(n)
        self.coords = newsies[-1]
        return newsies

    def get_potential(self, coords=None, atoms=None):
        """Handles potential calls

        :return: list of potentials
        :rtype:
        """
        if coords is None:
            coords = self.walkers.coords
        coords = np.asarray(coords)
        if atoms is None:
            atoms = self.coords.atoms

        if len(coords.shape) > 3:
            return np.array([ self.potential(atoms, crds) for crds in coords ])
        else:
            return self.potential(atoms, coords)

    @abstractmethod
    def branch(self):
        """Handles branching in DMC"""
        pass

    @abstractmethod
    def update_weights(self, potentials, v_refs):
        """Handles weights in DMC"""
        pass

    @abstractmethod
    def weight_descendants(self):
        """Applies descendant weighting to the history saved in self.stack"""
        pass

    def propagate(self, n=1):
        """Propagates the walker set n steps"""

        new_coords = self.diffuse(n)
        pots = self.get_potential(new_coords)
        v_refs = np.average(pots, axis=1)
        self.update_weights(pots, v_refs) # will write the weights to the WalkerSet at the end
        self.energies.append(v_refs) # we'll flatten this out at the end with itertools.chain
        self.snapshot()

    def snapshot(self, file=None):
        """Saves a snapshot of the simulation to file

        :param file:
        :type file:
        :return:
        :rtype:
        """
        import pickle

        if file is None:
            file = self.log
        if file is not None:
            with open(file, "w+") as dump:
                pickle.dump(self, dump)

