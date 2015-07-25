#!/usr/bin/env python

"""
Iteratively solve leaf temp, ci, gs, An and transpiration following Maetra
looping logic


That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.07.2015)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np
import os
import math

from farq import FarquharC3
from stomatal_conductance_models import StomtalConductance
from leaf_energy_balance import LeafEnergyBalance

class CoupledModel(object):
    """Iteratively solve leaf temp, ci, gs and An."""
    def __init__(self, g0, g1, D0, Vcmax25, Jmax25, Rd25, Eaj, Eav, deltaSj,
                 deltaSv, Hdv, Hdj, Q10, Ca, iter_max=100):

        # set params
        self.g0 = g0
        self.g1 = g1
        self.D0 = D0
        self.Vcmax25 = Vcmax25
        self.Jmax25 = Jmax25
        self.Rd25 = Rd25
        self.Eaj = Eaj
        self.Eav = Eav
        self.deltaSj = deltaSj
        self.deltaSv = deltaSv
        self.Hdv = Hdv
        self.Hdj = Hdj
        self.Q10 = Q10
        self.Ca = Ca
        self.leaf_width = leaf_width

        # Bit of hack to get around considering seperate sunlit/shaded leaves
        self.leaf_absorptance = leaf_absorptance
        self.iter_max = iter_max

        # Constants
        # Ratio of Gbh:Gbc
        self.GBHGBC = 1.32
        self.deg2kelvin = 273.15


    def main(self, tair, par, vpd, wind, pressure):

        F = FarquharC3(peaked_Jmax=True, peaked_Vcmax=False, model_Q10=True)
        S = StomtalConductance(g0=self.g0, g1=self.g1, D0=self.D0)
        L = LeafEnergyBalance(self.leaf_width, self.leaf_absorptance)

        # set initialise values
        dleaf = vpd
        Cs = Ca
        Tleaf = tair
        Tleaf_K = Tleaf + self.deg2kelvin

        print "Start: %.3f %.3f %.3f" % (Cs, Tleaf, dleaf)
        print

        iter = 0
        while True:

            (An, Acn,
             Ajn) = F.calc_photosynthesis(Ci=Cs, Tleaf=Tleaf_K, Par=par,
                                          Jmax25=Jmax25, Vcmax25=Vcmax25,
                                          Q10=Q10, Eaj=Eaj, Eav=Eav,
                                          deltaSj=deltaSj, deltaSv=deltaSv,
                                          Rd25=Rd25, Hdv=Hdv, Hdj=Hdj)
            gs = S.leuning(dleaf, An, Cs)

            (new_tleaf, et, gbH, gv) = L.calc_leaf_temp(Tleaf, tair, gs, par,
                                                        dleaf, pressure, wind)

            # update Cs and VPD
            gbc = gbH / self.GBHGBC
            Cs = Ca - An / gbc
            dleaf = et * pressure / gv

            print "%.3f %.3f %.3f %.3f %.3f" %  (Cs, Tleaf, dleaf, An, gs)

            if math.fabs(Tleaf - new_tleaf) < 0.02:
                break

            if iter > self.iter_max:
                raise Exception('No convergence: %d' % (iter))

            Tleaf = new_tleaf
            Tleaf_K = Tleaf + self.deg2kelvin
            iter += 1

        # Now recalculate new An and gs based on resolved vpd, ci, tleaf
        (An, Acn, Ajn) = F.calc_photosynthesis(Ci=Cs, Tleaf=Tleaf_K, Par=par,
                                               Jmax25=Jmax25, Vcmax25=Vcmax25,
                                               Q10=Q10, Eaj=Eaj,
                                               Eav=Eav, deltaSj=deltaSj,
                                               deltaSv=deltaSv, Rd25=Rd25,
                                               Hdv=Hdv, Hdj=Hdj)
        gs = S.leuning(dleaf, An, Cs)

        print
        print "End: %.3f %.3f %.3f %.3f %.3f" % (Cs, Tleaf, dleaf, An, gs)




if __name__ == '__main__':

    # Parameters

    # gs stuff
    g0 = 0.01
    g1 = 9.0
    D0 = 1.5 # kpa

    # A stuff
    Vcmax25 = 30.0
    Jmax25 = Vcmax25 * 2.0
    Rd25 = 2.0
    Eaj = 30000.0
    Eav = 60000.0
    deltaSj = 650.0
    deltaSv = 650.0
    Hdv = 200000.0
    Hdj = 200000.0
    Q10 = 2.0

    # Misc stuff
    leaf_width = 0.02

    # hack to get around doing seperate sunlit/shaded leaves
    leaf_absorptance = 0.5 # leaf absorptance of solar radiation [0,1]


    # variables though obviously fixed here.
    par = 1500.0
    tair = 20.0
    vpd = 2.0
    wind = 2.5
    pressure = 101.0
    Ca = 400.0

    C = CoupledModel(g0, g1, D0, Vcmax25, Jmax25, Rd25, Eaj, Eav, deltaSj,
                   deltaSv, Hdv, Hdj, Q10, Ca)
    C.main(tair, par, vpd, wind, pressure)
