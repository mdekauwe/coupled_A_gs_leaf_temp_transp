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
from solve_coupled_An_gs_leaf_temp_transpiration import CoupledModel
from utils import vpd_to_rh, get_dewpoint, calc_esat

def get_values(rh, Ca, tair, par, pressure, C):
    kpa_2_pa = 1000.
    pa_2_kpa = 1.0 / kpa_2_pa

    esat = calc_esat(tair, pressure)
    ea = rh / 100. * esat
    vpd = (esat - ea) * pa_2_kpa
    print rh, vpd
    gs_store = []
    et_store = []
    An_store = []
    tair_store = []
    for i,ta in enumerate(tair):

        #Td = get_dewpoint(ta, rh)
        #if Td > 0.0:
        (An, gsw, et, LE) = C.main(ta, par, vpd[i], wind, pressure, Ca)
        gs_store.append(gsw) # mol H20 m-2 s-1
        et_store.append(et*18*0.001*86400.) # mm d-1
        #et_store.append(LE)
        An_store.append(An) #umol m-2 s-1
        #An_store.append(An*12.*0.000001*86400.) # g C m-2 d-1
        tair_store.append(ta)
    return gs_store, et_store, An_store, tair_store

if __name__ == '__main__':

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.1)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['font.size'] = 10
    plt.rcParams['legend.fontsize'] = 8
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10

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



    # A stuff
    Vcmax25 = 40.0          # ENF CABLE value
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
    Ca2 = 900.

    g0 = 0.01
    g1 = 9.0
    D0 = 1.5
    gamma = 0.0
    CL = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, gs_model="leuning")

    g0 = 0.0
    g1 = 2.35
    CM = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, gs_model="medlyn")

    tair = np.linspace(0.1, 40, 50)



    #
    ## LEUNING
    #
    rh = 90.
    (gs_amb, et_amb, an_amb, tair_2plot) = get_values(rh, Ca1, tair, par, pressure, CL)
    (gs_ele, et_ele, an_ele, tair_2plot) = get_values(rh, Ca2, tair, par, pressure, CL)

    ax1.plot(tair_2plot, et_amb, "r-", label="LEU: %d (ppm)" % (int(Ca1)))
    ax1.plot(tair_2plot, et_ele, "r--", label="LEU: %d (ppm)" % (int(Ca2)))
    ax4.plot(tair_2plot, an_amb, "r-")
    ax4.plot(tair_2plot, an_ele, "r--")
    ax7.plot(tair_2plot, gs_amb, "r-")
    ax7.plot(tair_2plot, gs_ele, "r--")

    rh = 50.0
    (gs_amb, et_amb, an_amb, tair_2plot) = get_values(rh, Ca1, tair, par, pressure, CL)
    (gs_ele, et_ele, an_ele, tair_2plot) = get_values(rh, Ca2, tair, par, pressure, CL)

    ax2.plot(tair_2plot, et_amb, "r-")
    ax2.plot(tair_2plot, et_ele, "r--")
    ax5.plot(tair_2plot, an_amb, "r-")
    ax5.plot(tair_2plot, an_ele, "r--")
    ax8.plot(tair_2plot, gs_amb, "r-")
    ax8.plot(tair_2plot, gs_ele, "r--")

    rh = 10.0
    (gs_amb, et_amb, an_amb, tair_2plot) = get_values(rh, Ca1, tair, par, pressure, CL)
    (gs_ele, et_ele, an_ele, tair_2plot) = get_values(rh, Ca2, tair, par, pressure, CL)

    ax3.plot(tair_2plot, et_amb, "r-")
    ax3.plot(tair_2plot, et_ele, "r--")
    ax6.plot(tair_2plot, an_amb, "r-")
    ax6.plot(tair_2plot, an_ele, "r--")
    ax9.plot(tair_2plot, gs_amb, "r-")
    ax9.plot(tair_2plot, gs_ele, "r--")


    #
    ## MEDLYN
    #


    rh = 90.0
    (gs_amb, et_amb, an_amb, tair_2plot) = get_values(rh, Ca1, tair, par, pressure, CM)
    (gs_ele, et_ele, an_ele, tair_2plot) = get_values(rh, Ca2, tair, par, pressure, CM)

    ax1.plot(tair_2plot, et_amb, "g-", label="MED: %d (ppm)" % (int(Ca1)))
    ax1.plot(tair_2plot, et_ele, "g--", label="MED: %d (ppm)" % (int(Ca2)))
    ax4.plot(tair_2plot, an_amb, "g-")
    ax4.plot(tair_2plot, an_ele, "g--")
    ax7.plot(tair_2plot, gs_amb, "g-")
    ax7.plot(tair_2plot, gs_ele, "g--")

    rh = 50.0
    (gs_amb, et_amb, an_amb, tair_2plot) = get_values(rh, Ca1, tair, par, pressure, CM)
    (gs_ele, et_ele, an_ele, tair_2plot) = get_values(rh, Ca2, tair, par, pressure, CM)

    ax2.plot(tair_2plot, et_amb, "g-")
    ax2.plot(tair_2plot, et_ele, "g--")
    ax5.plot(tair_2plot, an_amb, "g-")
    ax5.plot(tair_2plot, an_ele, "g--")
    ax8.plot(tair_2plot, gs_amb, "g-")
    ax8.plot(tair_2plot, gs_ele, "g--")

    rh = 10.0
    (gs_amb, et_amb, an_amb, tair_2plot) = get_values(rh, Ca1, tair, par, pressure, CM)
    (gs_ele, et_ele, an_ele, tair_2plot) = get_values(rh, Ca2, tair, par, pressure, CM)

    ax3.plot(tair_2plot, et_amb, "g-")
    ax3.plot(tair_2plot, et_ele, "g--")
    ax6.plot(tair_2plot, an_amb, "g-")
    ax6.plot(tair_2plot, an_ele, "g--")
    ax9.plot(tair_2plot, gs_amb, "g-")
    ax9.plot(tair_2plot, gs_ele, "g--")

    ax8.set_xlabel("Air temperature ($^{\circ}$C)")
    ax1.set_ylabel("$E$ (mm d$^{-1}$)")
    ax4.set_ylabel("$A_{\mathrm{n}}$ ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)")
    ax7.set_ylabel("$g_{\mathrm{s}}$ (mol m$^{-2}$ s$^{-1}$)")

    ax1.get_yaxis().set_label_coords(-0.15,0.5)
    ax4.get_yaxis().set_label_coords(-0.15,0.5)
    ax7.get_yaxis().set_label_coords(-0.15,0.5)

    ax1.legend(numpoints=1, loc="best", frameon=False)

    ax1.set_ylim(0,4)
    ax2.set_ylim(0,4)
    ax3.set_ylim(0,4)
    ax4.set_ylim(0,15)
    ax5.set_ylim(0,15)
    ax6.set_ylim(0,15)
    ax7.set_ylim(0,0.18)
    ax8.set_ylim(0,0.18)
    ax9.set_ylim(0,0.18)


    ax1.set_xlim(0,40)
    ax2.set_xlim(0,40)
    ax3.set_xlim(0,40)
    ax4.set_xlim(0,40)
    ax5.set_xlim(0,40)
    ax6.set_xlim(0,40)
    ax7.set_xlim(0,40)
    ax8.set_xlim(0,40)
    ax9.set_xlim(0,40)

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

    ax1.set_title("RH $=$ 90 %")
    ax2.set_title("RH $=$ 50 %")
    ax3.set_title("RH $=$ 10 %")

    fig.savefig("/Users/mdekauwe/Desktop/Fig_S01.pdf", bbox_inches='tight',
                pad_inches=0.1)
    plt.show()
