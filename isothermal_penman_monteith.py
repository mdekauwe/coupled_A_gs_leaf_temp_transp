#!/usr/bin/env python

"""
Isothermal Penman-Monteith

Reference:
==========
* Leuning et al. (1995) Leaf nitrogen, photosynthesis, conductance and
  transpiration: scaling from leaves to canopies

"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.07.2015)"
__email__ = "mdekauwe@gmail.com"


import math
import sys


class PenmanMonteith(object):

    """docstring for PenmanMonteith"""
    def __init__(self, leaf_width, leaf_absorptance):

        self.kpa_2_pa = 1000.
        self.sigma = 5.6704E-08  # stefan boltzmann constant, (w m-2 k-4)
        self.emissivity_leaf = 0.99   # emissivity of leaf (-)
        self.cp = 1010.0         # specific heat of dry air (j kg-1 k-1)
        self.h2olv0 = 2.501e6    # latent heat H2O (J kg-1)
        self.h2omw = 18E-3       # mol mass H20 (kg mol-1)
        self.air_mass = 29.0E-3     # mol mass air (kg mol-1)
        self.umol_to_j = 4.57    # conversion from J to umol quanta
        self.dheat = 21.5E-6     # molecular diffusivity for heat (m2 s-1)
        self.DEG_TO_KELVIN = 273.15
        self.RGAS = 8.314        # universal gas constant (mol-1 K-1)
        self.leaf_absorptance = leaf_absorptance
        self.leaf_width = leaf_width # (m)

    def calc_et(self, tleaf, tair, gs, vpd, pressure, wind, par, gh, gv,
                rnet):

        # latent heat of water vapour at air temperature (j mol-1)
        lambda_et = (self.h2olv0 - 2.365e3 * tair) * self.h2omw

        # curve relating sat water vapour pressure to temperature (pa K-1)
        # kelvin conversion in func
        slope = ((self.calc_esat(tair + 0.1, pressure) -
                  self.calc_esat(tair, pressure)) / 0.1)

        # psychrometric constant
        gamma = self.cp * self.air_mass * pressure * 1000.0 / lambda_et

        # Y cancels in eqn 10
        arg1 = (slope * rnet + (vpd * self.kpa_2_pa) * gh * self.cp *
                self.air_mass)
        arg2 = slope + gamma * gh / gv
        et = arg1 / arg2

        # latent heat loss
        LE_et = et

        # et units = mol m-2 s-1,
        # multiply by 18 (grams)* 0.001 (grams to kg) * 86400.
        # to get to kg m2 d-1 or mm d-1
        return et / lambda_et, LE_et

    def calc_conductances(self, tair_k, tleaf, tair, pressure, wind, gs,
                          cmolar):
        """
        Both forced and free convection contribute to exchange of heat and mass
        through leaf boundary layers at the wind speeds typically encountered
        within plant canopies (<0-5ms~'). It is particularly imponant to
        includethe contribution of buoyancy forces to the boundary
        conductance for sunlit leaves deep within the canopy where wind speeds
        are low, for without this mechanism computed leaf temperatures become
        excessively high.

        Leuning 1995, appendix E
        Medlyn et al. 2007 appendix, for need for cmolar
        """
        # Ratio of Gbw:Gbh
        GBVGBH = 1.075

        # Ratio of Gsw:Gsc
        GSVGSC = 1.57

        # radiation conductance (mol m-2 s-1)
        grn = ((4.0 * self.sigma * tair_k**3 * self.emissivity_leaf) /
               (self.cp * self.air_mass))

        # boundary layer conductance for 1 side of leaf from forced convection
        # (mol m-2 s-1)
        gbHw = 0.003 * math.sqrt(wind / self.leaf_width) * cmolar

        # grashof number
        grashof_num = 1.6e8 * math.fabs(tleaf - tair) * self.leaf_width**3

        # boundary layer conductance for free convection
        # (mol m-2 s-1)
        gbHf = 0.5 * self.dheat * grashof_num**0.25 / self.leaf_width * cmolar

        # total boundary layer conductance to heat for one side of the leaf
        gbH = gbHw + gbHf

        # ... for hypostomatous leaves only gbH should be doubled and the
        # single-sided value used for gbw

        # total leaf conductance to heat (mol m-2 s-1), two sided see above.
        gh = 2.0 * (gbH + grn)

        gbv = GBVGBH * gbH
        gsv = GSVGSC * gs

        # total leaf conductance to water vapour (mol m-2 s-1)
        gv = (gbv * gsv) / (gbv + gsv)

        return (grn, gh, gbH, gv)

    def calc_rnet(self, pressure, par, tair, tair_k, tleaf_k, vpd):

        umol_m2_s_to_W_m2 = 2.0 / self.umol_to_j
        par *= umol_m2_s_to_W_m2

        # atmospheric water vapour pressure (Pa)
        ea = max(0.0, self.calc_esat(tair, pressure) - (vpd * self.kpa_2_pa))

        # eqn D4
        emissivity_atm = 0.642 * (ea / tair_k)**(1.0 / 7.0)

        rlw_down = emissivity_atm * self.sigma * tair_k**4
        rlw_up = self.emissivity_leaf * self.sigma * tleaf_k**4
        isothermal_net_lw = rlw_up - rlw_down

        # isothermal net radiation (W m-2)
        return (self.leaf_absorptance * par - isothermal_net_lw)

    def calc_slope_of_saturation_vapour_pressure_curve(self, tair):
        """ Eqn 13 from FAO paper, Allen et al. 1998.

        Returns:
        --------
        slope : float
            slope of saturation vapour pressure curve [kPa degC-1]

        """
        t = tair + 237.3
        arg1 = 4098.0 * (0.6108 * math.exp((17.27 * tair) / t))
        arg2 = t**2
        return (arg1 / arg2) * self.kpa_2_pa

    def calc_esat(self, temp, pressure):
        """
        Saturation vapor pressure (Pa K-1)

        Values of saturation vapour pressure from the Tetens formula are
        within 1 Pa of the exact values.

        Taken from Stull 2000 Meteorology for Scientist and Engineers, but see
        also Jones 1992 p 110 (note error in a - wrong units)
        """
        Tk = temp + self.DEG_TO_KELVIN
        e0 = 0.611 * self.kpa_2_pa
        b = 17.2694
        T1 = 273.16 # kelvin
        T2 = 35.86  # kelvin

        return e0 * math.exp(b * (Tk - T1) / (Tk - T2))

    def main(self, tleaf, tair, gs, vpd, pressure, wind, par):
        tleaf_k = tleaf + DEG_TO_KELVIN
        tair_k = tair + DEG_TO_KELVIN

        # density of dry air
        air_density = pressure * 1000.0 / (287.058 * tair_k)
        cmolar = pressure * 1000.0 / (RGAS * tair_k)

        rnet_iso = P.calc_rnet(pressure, par, tair, tair_k, tleaf_k, vpd)

        (grn, gh, gbH, gv) = P.calc_conductances(tair_k, tleaf, tair, pressure,
                                                wind, gs, cmolar)
        (et, lambda_et) = P.calc_et(tleaf, tair, gs, vpd, pressure, wind, par,
                                    gh, gv, rnet_iso)
        return (et, lambda_et)

if __name__ == '__main__':

    par = 1000.0
    tleaf = 21.5
    tair = 20.0
    gs = 0.15
    vpd = 2.0
    pressure = 101.0 # Pa
    wind = 2.0
    leaf_width = 0.02
    leaf_absorptance = 0.5 # leaf absorptance of solar radiation [0,1]
    DEG_TO_KELVIN = 273.15
    RGAS = 8.314

    P = PenmanMonteith(leaf_width, leaf_absorptance)
    (et, lambda_et) = P.main(tleaf, tair, gs, vpd, pressure, wind, par)

    print et, lambda_et
