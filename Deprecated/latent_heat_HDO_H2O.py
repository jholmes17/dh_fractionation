################################################################################
# find_LDO_expression.py
# TYPE: Supporting
# WHICH: Equilibrium and perturbation experiments
# DESCRIPTION: This script uses data from Jancso 1974
# (https://pubs.acs.org/doi/pdf/10.1021/cr60292a004) and Jacobson's Fundamentals
# of Atmospheric Modeling to get an analytical expression for L_HDO at
# temperatures in the Mars range of interest (70K to 350K).
#
# How it does it:
#     1. Used Jancso 1974, pg. 734, to get data on the difference of latent heat
#        of evaporation (L_e) for HDO and H2O from 0 to 100°C.
#     2. Used equation 2.54 in Jacobson to get analytical values for L_e of H2O
#        in the same range (0-100°C)
#     3. Combine the L_e of H2O with data from Jancso 1974 to get L_HDO.
#     4. Noting that L_e expression is roughly linear, fit a line to the L_HDO
#        values and extrapolated to Mars temperatures of interest.
#
# Eryn Cangi
# Originally written 2018, updated to be scientifically respectable 1 July 2019
# Currently tested for Python: 3.6
# NOTE 2 SEPT 2019:  THIS IS DEPRECATED. This calculation is not necessary
################################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def Lguess_LINEAR(t, M, a, B):
    return M*(t-a) + B


def L_evap(Tc):
    L_JpKg = 2.501e6 - 2370*Tc
    L_calpmol = L_JpKg * (1/4.184) * (1/1000) * (18.015/1)
    return L_calpmol


def KtoC(T):
    return T - 273.15


jancso_data = np.array([219, 182, 151, 125, 102])

temps = np.asarray(range(0, 125, 25))
L_H2O = np.array([L_evap(t) for t in temps])
L_HDO = L_H2O + jancso_data

popt, pcov = curve_fit(Lguess_LINEAR, temps, L_HDO, p0=[-11, 0, 11000])
a = popt[0]
b = popt[1]
c = popt[2]

marsrange = np.arange(KtoC(70), KtoC(375))

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Louis George Caf?'
plt.rcParams['font.monospace'] = 'FreeMono'
plt.rcParams['font.size'] = 19
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22

fig = plt.figure(figsize=(8, 6))
ax1 = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.3)

# Celsius axis
ax1.grid()
ax1.set_xticks([-200, -100, 0, 100])
ax1.scatter(temps, L_HDO, label="HDO, Jancso+1974")
ax1.plot(marsrange, [L_evap(t) for t in marsrange], label="H$_2$O")
ax1.plot(marsrange, Lguess_LINEAR(marsrange, a, b, c), label="linear fit",
         color="xkcd:forest green")
ax1.set_ylabel("L$_e$ (cal/mol)")
ax1.set_xlabel("Temperature (C)")
ax1.legend()

# Kelvin axis
ax2 = ax1.twiny()
ax2.spines["bottom"].set_position(("axes", -0.3))
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")
ax2.set_frame_on(True)
ax2.patch.set_visible(False)
for side in ["top", "left", "right"]:
    ax2.spines[side].set_visible(False)
ax2.set_xticks([-200, -100, 0, 100])
ax2.set_xticklabels([round(tic + 273.15) for tic in [-200, -100, 0, 100]])
ax2.set_xlabel("Temperature (K)")
ax2.margins(ax1.margins()[0], ax1.margins()[1])
plt.title("H$_2$O and HDO latent heat of evaporation")
plt.show()
