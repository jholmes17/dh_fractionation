################################################################################
# h2oloss_family.jl
# TYPE: MAIN (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: make plot of loss of H2O in a family of solutions with
# various D/H ratios and fractionation factor. NOTE: This file assumes the 
# current water inventory is 25 m GEL, per Villanueva 2015; one could assume
# 20 or 30 and get different answers as these are all relatively in agreement 
# with Villanueva 2015 and other geological estimates (e.g. Carr 2019 "Formation
# and fate of a frozen Hesperian ocean".
#
# Eryn Cangi
# 2019
# Currently tested for Python: 3.7
################################################################################

import numpy as np
import matplotlib.pyplot as plt


def h2olost(cur_h2o, DHnow=5.5, DH0=1, f=0.001):
    return cur_h2o * ((DHnow/DH0)**(1/(1-f)) - 1)


def plot_water_loss_results():
    current_inventory = 25
    farray = np.linspace(0.000005, 0.5, 100)
    f_thermal = [0.0012, 0.0012, 0.016, 0.000043, 0.0022, 0.00055, 0.0013,
                 0.0011, 0.0011]
    f_nonthermal = [0.038, 0.032, 0.065, 0.038, 0.068, 0.023, 0.037,
                    0.037, 0.039]
    loss_rd = [h2olost(current_inventory, 5.5, 1.275, f) for f in farray]
    loss_thermal = [h2olost(current_inventory, 5.5, 1.275, f) for f in f_thermal]
    loss_nonthermal = [h2olost(current_inventory, 5.5, 1.275, f) for f in f_nonthermal]

    # make plot
    fig = plt.figure(figsize=(10, 7))

    plt.style.use('default')
    plt.rc('text', usetex=False)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Louis George Caf?'
    plt.rcParams['font.monospace'] = 'FreeMono'
    plt.rcParams['font.size'] = 26
    plt.rcParams['axes.labelsize'] = 28
    plt.rcParams['xtick.labelsize'] = 26
    plt.rcParams['ytick.labelsize'] = 26

    ax = plt.gca()
    ax.set_facecolor("#ededed")
    ax.grid(zorder=0, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]:
        ax.spines[side].set_visible(False)

    plt.xlabel(r"Fractionation factor $f$")
    plt.ylabel("GEL water lost (m) since 4.5 Ga")
    plt.title("Water loss (Rayleigh distillation)")

    # plot things
    deemph = "#555555"
    # ax.semilogx(farray, l1, zorder=2, color="xkcd:vivid purple")
    ax.semilogx(farray, loss_rd, zorder=2, color="cornflowerblue")
    # ax.semilogx(farray, l3, zorder=2, color="orange")
    # ax.semilogx(farray, l4, zorder=2, color="red")
    ax.scatter(f_thermal, loss_thermal, zorder=3, color="black", marker="x", s=150,
               label="This work, thermal escape only")
    ax.scatter(f_nonthermal, loss_nonthermal, zorder=3, color="black", marker="d", s=150,
               label="This work, thermal+nonthermal escape")
    ax.scatter(0.32, 3.4, zorder=2, marker="v", s=100, color=deemph) # Yung
    ax.scatter(0.016, 65, zorder=2, marker="*", s=100, color=deemph) # Kras 2000 Model 1
    ax.scatter(0.135, 160, zorder=2,marker="^", s=100, color=deemph) # Kras 2000 Model 2

    # Plot annotations
    # plt.text(0.0000135, 40, r"(D/H)$_{M}$ = $3 \times$(D/H)$_{\oplus}$", color="xkcd:vivid purple", fontsize=22))
    plt.text(0.0000135, 90, r"(D/H)$_{M}$ = $5.5 \times$(D/H)$_{E}$", color="cornflowerblue", fontsize=22)
    # plt.text(0.0000135, 140, r"(D/H)$_{M}$ = $8 \times$(D/H)$_{\oplus}$", color="orange", fontsize=22))
    # plt.text(0.0000135, 180, r"(D/H)$_{M}$ = $10 \times$(D/H)$_{\oplus}$", color="red", fontsize=22))
    plt.text(0.1, 150, "K2000, Model 1", ha="right", color=deemph, fontsize=22)
    plt.text(0.02, 60, "K2000, Model 2", color=deemph, fontsize=22)
    plt.text(0.03, 0, "Yung 1988", color=deemph, fontsize=22)
    plt.legend(fontsize=24)
    plt.ylim(0, 240)
    plt.savefig("../Results/ALL STUDY PLOTS/h2oloss.png", bbox_inches="tight")
    print("Water loss values:")
    print(loss_nonthermal)
    print(loss_thermal)


def plot_simple_waterloss():
    current_inventory = 25
    farray = np.linspace(0.000005, 0.5, 100)

    l1 = [h2olost(current_inventory, 5.5, 1.275, f) for f in farray]

    # make plot
    fig = plt.figure(figsize=(10, 7))

    plt.style.use('default')
    plt.rc('text', usetex=False)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Louis George Caf?'
    plt.rcParams['font.monospace'] = 'FreeMono'
    plt.rcParams['font.size'] = 26
    plt.rcParams['axes.labelsize'] = 28
    plt.rcParams['xtick.labelsize'] = 26
    plt.rcParams['ytick.labelsize'] = 26

    ax = plt.gca()
    ax.set_facecolor("#ededed")
    ax.grid(zorder=0, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]:
        ax.spines[side].set_visible(False)

    plt.xlabel(r"Fractionation factor $f$")
    plt.ylabel("Watr lost (m GEL) in 4.5 Ga")

    # plot things
    ax.semilogx(farray, l1, zorder=2, color="cornflowerblue")
    # plt.legend(fontsize=24)
    plt.text(0.000005, 50, "Water loss constant \n"+r" with $f$ when $f\ll0.05$")
    plt.text(0.0006, 200, "Water loss depends \nstrongly on " + r"$f$ for $f\gtrapprox 0.05$")
    plt.ylim(0, 240)
    plt.savefig("../Results/ALL STUDY PLOTS/h2oloss_simple.png", bbox_inches="tight")


plot_water_loss_results()
plot_simple_waterloss()