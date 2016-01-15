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

def get_values(vpd_vals, Ca, tair, par, pressure, C):

    gs_store = []
    et_store = []
    An_store = []
    for i,vpd in enumerate(vpd_vals):

        #Td = get_dewpoint(ta, rh)
        #if Td > 0.0:
        (An, gsw, et) = C.main(tair, par, vpd, wind, pressure, Ca)
        gs_store.append(gsw)
        et_store.append(et*18*0.001*86400.)
        An_store.append(An*12.*0.000001*86400.)
    return gs_store, et_store, An_store

if __name__ == '__main__':

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.1)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
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

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)



    # Parameters



    # A stuff
    Vcmax25 = 35.0          # ENF CABLE value
    Jmax25 = Vcmax25 * 2.0
    Rd25 = 0.54
    Eaj = 37870.0
    Eav = 60880.0
    deltaSj = 640.0
    deltaSv = 640.0
    Hdv = 200000.0
    Hdj = 200000.0
    Q10 = 2.0
    alpha = 0.385
    # Misc stuff
    leaf_width = 0.01
    SW_abs = 0.5 # absorptance to short_wave rad [0,1], typically 0.4-0.6


    # variables though obviously fixed here.
    par = 1500.0
    wind = 2.5
    pressure = 101325.0
    tair = 25.0
    Ca = 380.
    vpd = np.linspace(0.14, 9, 100)

    g0 = 0.0
    g1 = 9.0
    D0 = 1.5
    gamma = 0.0
    CL = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, alpha=alpha, gs_model="leuning")

    g0 = 0.0
    g1 = 3.0
    CM = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, alpha=alpha, gs_model="medlyn")


    (gs, et, an) = get_values(vpd, Ca, tair, par, pressure, CM)
    ax1.plot(vpd, et, "b-", label="MED")
    ax3.plot(vpd, gs, "b-", label="MED")

    (gs, et, an) = get_values(vpd, Ca, tair, par, pressure, CL)
    ax1.plot(vpd, et, "r-", label="LEU")
    ax3.plot(vpd, gs, "r-", label="LEU")



    g0 = 0.01
    g1 = 9.0
    D0 = 1.5
    gamma = 0.0
    CLX = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, alpha=alpha, gs_model="leuning")

    g0 = 0.01
    g1 = 3.0
    CMX = CoupledModel(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25,
                     Eaj, Eav,deltaSj, deltaSv, Hdv, Hdj, Q10, leaf_width,
                     SW_abs, alpha=alpha, gs_model="medlyn")

    (gs, et, an) = get_values(vpd, Ca, tair, par, pressure, CLX)
    ax2.plot(vpd, et, "r-", label="LEU")
    ax4.plot(vpd, gs, "r-", label="LEU")

    (gs, et, an) = get_values(vpd, Ca, tair, par, pressure, CMX)
    ax2.plot(vpd, et, "b-", label="MED")
    ax4.plot(vpd, gs, "b-", label="MED")


    ax1.legend(numpoints=1, loc="best", frameon=False)

    ax1.set_ylabel("$E$ (mm d$^{-1}$)")
    ax3.set_ylabel("$g_{\mathrm{s}}$ (mol m$^{-2}$ s$^{-1}$)")
    ax3.set_xlabel("$D$ (kPa)", position=(1.0, 0.5))

    ax1.set_ylim(0,6)
    ax2.set_ylim(0,6)
    ax3.set_ylim(0,0.18)
    ax4.set_ylim(0,0.18)


    ax1.locator_params(nbins=6)
    ax2.locator_params(nbins=6)
    ax3.locator_params(nbins=6)
    ax4.locator_params(nbins=6)

    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    plt.show()


    #fig.savefig("/Users/%s/Desktop/Fig_S01.pdf" % (os.getlogin()),
    #            bbox_inches='tight', pad_inches=0.1)
