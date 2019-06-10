from .ContinuousWeightingDMC import ContinuousWeightingDMC

class ImportanceSamplingDMC(ContinuousWeightingDMC):

    def get_potential(self, coords=None, atoms=None):
        # this needs an update here
        raise NotImplemented