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
import pandas as pd


def h2olost(cur_h2o, DHnow=5.5, DH0=1, f=0.001):
    return cur_h2o * ((DHnow/DH0)**(1/(1-f)) - 1)


current_inventory = 25
farray = np.linspace(0.000005, 0.35, 100)
farray2 = [0.00036, 0.00089, 0.00064, 0.00016, 0.00078, 1.16e-5, 0.0036,
           0.00037, 0.00036]
l1 = [h2olost(current_inventory, 3, 1.275, f) for f in farray]
l2 = [h2olost(current_inventory, 5.5, 1.275, f) for f in farray]
l3 = [h2olost(current_inventory, 8, 1.275, f) for f in farray]
l4 = [h2olost(current_inventory, 10, 1.275, f) for f in farray]

model_output = [h2olost(current_inventory, 5.5, 1.275, f2) for f2 in farray2]

# make plot
fig = plt.figure(figsize=(11, 8))

plt.style.use('default')
plt.rc('text', usetex=False)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Louis George Caf?'
plt.rcParams['font.monospace'] = 'FreeMono'
plt.rcParams['font.size'] = 22
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22

ax = plt.gca()
ax.set_facecolor("#ededed")
ax.grid(zorder=0, color="white", which="major")
for side in ["top", "bottom", "left", "right"]:
    ax.spines[side].set_visible(False)

plt.xlabel(r"Fractionation factor $f$")
plt.ylabel("GEL water lost (m) since 4.5 Ga")
#plt.title("Water loss vs. fractionation factor", fontsize=36)

# plot things
ax.semilogx(farray, l1, zorder=2, color="purple")
ax.semilogx(farray, l2, zorder=2, color="cornflowerblue")
ax.semilogx(farray, l3, zorder=2, color="orange")
ax.semilogx(farray, l4, zorder=2, color="red")
ax.scatter(farray2, model_output, zorder=3, color="black", marker="x", s=100,
           label="This work")
ax.scatter(0.32, 3.4, zorder=2, color="black", marker="v",  s=100) # Yung
ax.scatter(0.016, 65, zorder=2, color="black", marker="*", s=100) # Kras 2000 Model 1
ax.scatter(0.135, 160, zorder=2, color="black", marker="^", s=100) # Kras 2000 Model 2


# Plot annotations
plt.text(0.0000135, 40, r"D/H = $3 \times$ SMOW", color="purple")
plt.text(0.0000135, 90, r"D/H = $5.5 \times$ SMOW", color="cornflowerblue")
plt.text(0.0000135, 140, r"D/H = $8 \times$ SMOW", color="orange")
plt.text(0.0000135, 180, r"D/H = $10 \times$ SMOW", color="red")
plt.text(0.000005, 250, r"For very small $f$," + " water \nloss only depends on \ncurrent D/H ratio and \ncurrent water inventory, \nnot f.")
plt.text(0.005, 400, r"When $f>$0.05, " + "\nwater loss depends \nmuch more strongly \non f.")
plt.text(0.1, 150, "K2000, Model 1", ha="right")
plt.text(0.02, 60, "K2000, Model 2")
plt.text(0.035, 0, "Yung 1988")
plt.legend(fontsize=18)
plt.savefig("h2oloss_family.png", bbox_inches="tight")

