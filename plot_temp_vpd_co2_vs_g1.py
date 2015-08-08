#!/usr/bin/env python

"""
Compare transpiration sensitivity to PFT difference in g1 vs. inc. temp/VPD

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.07.2015)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np
import os
import math
import matplotlib.pyplot as plt


from farq import FarquharC3
from leaf_energy_balance import LeafEnergyBalance
from solve_coupled_An_gs_leaf_temp_transpiration import CoupledModel


def get_values(vpd, Ca, tair, par, pressure, C):
    gs_store = []
    et_store = []
    An_store = []
    for ta in tair:
        (An, gs, et) = C.main(ta, par, vpd, wind, pressure, Ca)
        gs_store.append(gs)
        et_store.append(et*18*0.001*86400.)
        An_store.append(An*12.*0.000001*86400.)
    return gs_store, et_store, An_store

if __name__ == '__main__':

    fig = plt.figure(figsize=(12,9))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.1)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color'] = almost_black
    plt.rcParams['xtick.color'] = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color'] = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor'] = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    #colour_list = brewer2mpl.get_map('Accent', 'qualitative', 8).mpl_colors
    # CB palette  with grey:
    # from http://jfly.iam.u-tokyo.ac.jp/color/image/pallete.jpg
    colour_list = ["#CC79A7", "#E69F00", "#0072B2", "#009E73", "#F0E442",
                   "#56B4E9", "#D55E00", "#000000"]

    ax1 = fig.add_subplot(331)
    ax2 = fig.add_subplot(332)
    ax3 = fig.add_subplot(333)
    ax4 = fig.add_subplot(334)
    ax5 = fig.add_subplot(335)
    ax6 = fig.add_subplot(336)
    ax7 = fig.add_subplot(337)
    ax8 = fig.add_subplot(338)
    ax9 = fig.add_subplot(339)


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
    leaf_width = 0.01

    SW_abs = 0.5 # absorptance to short_wave rad [0,1], typically 0.4-0.6


    # variables though obviously fixed here.
    par = 1500.0

    wind = 2.5
    pressure = 101325.0

    Ca1 = 500.
    Ca2 = 800.

    tair = np.linspace(0, 40, 40)
    C = CoupledModel(g0, g1, D0, Vcmax25, Jmax25, Rd25, Eaj, Eav, deltaSj,
                     deltaSv, Hdv, Hdj, Q10, leaf_width, SW_abs,
                     gs_model="leuning")
    vpd = 1.0
    (gs_amb, et_amb, an_amb) = get_values(vpd, Ca1, tair, par, pressure, C)
    (gs_ele, et_ele, an_ele) = get_values(vpd, Ca2, tair, par, pressure, C)

    ax1.plot(tair, et_amb, "r-", label="LEU: %d (ppm)" % (int(Ca1)))
    ax1.plot(tair, et_ele, "r--", label="LEU: %d (ppm)" % (int(Ca2)))
    ax4.plot(tair, an_amb, "r-")
    ax4.plot(tair, an_ele, "r--")
    ax7.plot(tair, gs_amb, "r-")
    ax7.plot(tair, gs_ele, "r--")

    vpd = 3.0
    gs_amb, et_amb, an_amb = get_values(vpd, Ca1, tair, par, pressure, C)
    gs_ele, et_ele, an_ele = get_values(vpd, Ca2, tair, par, pressure, C)

    ax2.plot(tair, et_amb, "r-")
    ax2.plot(tair, et_ele, "r--")
    ax5.plot(tair, an_amb, "r-")
    ax5.plot(tair, an_ele, "r--")
    ax8.plot(tair, gs_amb, "r-")
    ax8.plot(tair, gs_ele, "r--")

    vpd = 5.0
    gs_amb, et_amb, an_amb = get_values(vpd, Ca1, tair, par, pressure, C)
    gs_ele, et_ele, an_ele = get_values(vpd, Ca2, tair, par, pressure, C)

    ax3.plot(tair, et_amb, "r-")
    ax3.plot(tair, et_ele, "r--")
    ax6.plot(tair, an_amb, "r-")
    ax6.plot(tair, an_ele, "r--")
    ax9.plot(tair, gs_amb, "r-")
    ax9.plot(tair, gs_ele, "r--")

    

    g0 = 0.01
    g1 = 2.35
    D0 = -999.9 # kpa
    C = CoupledModel(g0, g1, D0, Vcmax25, Jmax25, Rd25, Eaj, Eav, deltaSj,
                     deltaSv, Hdv, Hdj, Q10, leaf_width, SW_abs,
                     gs_model="medlyn")


    vpd = 1.0
    gs_amb, et_amb, an_amb = get_values(vpd, Ca1, tair, par, pressure, C)
    gs_ele, et_ele, an_ele = get_values(vpd, Ca2, tair, par, pressure, C)

    ax1.plot(tair, et_amb, "g-", label="MED: %d (ppm)" % (int(Ca1)))
    ax1.plot(tair, et_ele, "g--", label="MED: %d (ppm)" % (int(Ca2)))
    ax4.plot(tair, an_amb, "g-")
    ax4.plot(tair, an_ele, "g--")
    ax7.plot(tair, gs_amb, "g-")
    ax7.plot(tair, gs_ele, "g--")

    vpd = 3.0
    gs_amb, et_amb, an_amb = get_values(vpd, Ca1, tair, par, pressure, C)
    gs_ele, et_ele, an_ele = get_values(vpd, Ca2, tair, par, pressure, C)

    ax2.plot(tair, et_amb, "g-")
    ax2.plot(tair, et_ele, "g--")
    ax5.plot(tair, an_amb, "g-")
    ax5.plot(tair, an_ele, "g--")
    ax8.plot(tair, gs_amb, "g-")
    ax8.plot(tair, gs_ele, "g--")

    vpd = 5.0
    gs_amb, et_amb, an_amb = get_values(vpd, Ca1, tair, par, pressure, C)
    gs_ele, et_ele, an_ele = get_values(vpd, Ca2, tair, par, pressure, C)

    ax3.plot(tair, et_amb, "g-")
    ax3.plot(tair, et_ele, "g--")
    ax6.plot(tair, an_amb, "g-")
    ax6.plot(tair, an_ele, "g--")
    ax9.plot(tair, gs_amb, "g-")
    ax9.plot(tair, gs_ele, "g--")

    ax8.set_xlabel("Tair ($^{\circ}$C)")
    ax1.set_ylabel("$E$ (mm d$^{-1}$)")
    ax4.set_ylabel("$A_{\mathrm{n}}$ (g C m$^{-2}$ d$^{-1}$)")
    ax7.set_ylabel("$g_{\mathrm{s}}$ (mmol m$^{-2}$ s$^{-1}$)")
    ax1.legend(numpoints=1, loc="best", frameon=False)

    ax1.set_ylim(0,4)
    ax2.set_ylim(0,4)
    ax3.set_ylim(0,4)
    ax4.set_ylim(0,12)
    ax5.set_ylim(0,12)
    ax6.set_ylim(0,12)
    ax7.set_ylim(0,0.06)
    ax8.set_ylim(0,0.06)
    ax9.set_ylim(0,0.06)

    ax1.locator_params(nbins=4)
    ax2.locator_params(nbins=4)
    ax3.locator_params(nbins=4)
    ax4.locator_params(nbins=6)
    ax5.locator_params(nbins=6)
    ax6.locator_params(nbins=6)
    ax7.locator_params(nbins=6)
    ax8.locator_params(nbins=6)
    ax9.locator_params(nbins=6)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax6.get_xticklabels(), visible=False)

    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)
    plt.setp(ax6.get_yticklabels(), visible=False)
    plt.setp(ax8.get_yticklabels(), visible=False)
    plt.setp(ax9.get_yticklabels(), visible=False)

    ax1.set_title("VPD = 1.0 (kPa)")
    ax2.set_title("VPD = 3.0 (kPa)")
    ax3.set_title("VPD = 5.0 (kPa)")

    fig.savefig("/Users/mdekauwe/Desktop/temp_co2_vs_g1.pdf", bbox_inches='tight',
                pad_inches=0.1)
    plt.show()
