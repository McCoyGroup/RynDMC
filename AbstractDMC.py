from abc import *
from .WalkerSet import WalkerSet


##################################################################################################################
##
##                                                  AbstractDMC
##
##
class AbstractDMC(metaclass=ABCMeta):
    """A completely abstract implementation of """

    def __init__(self, walker_set):
        if not isinstance(walker_set, WalkerSet):
            walker_set=WalkerSet(*walker_set)
        self.walkers = walker_set

    @abstractmethod
    def branch(self):
        """Handles branching in DMC"""
        pass

    @abstractmethod
    def propagate(self, n=1):
        """Propagates the walker set one step"""
        pass
