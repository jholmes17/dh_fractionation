################################################################################
# f_barchart.jl
# TYPE: MAIN (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: plots results for f in equilibrated atmosphere simulations in a
# bar chart, by experiment. Done in Python because it doesn't need to call on 
# Julia data/output.
#
# Eryn Cangi
# July 2019
# Currently tested for Python: 3.7
################################################################################
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rc
rc('text', usetex=True)

experiments = pd.DataFrame([["Mean", 3.6e-4],
                            [r"$\overline{T}_{surf}-25$%", 8.9e-4],
                            [r"$\overline{T}_{surf}+25$%", 6.4e-4],
                            [r"$\overline{T}_{tropo}-25$%", 1.6e-4],
                            [r"$\overline{T}_{tropo}+25$%", 7.8e-4],
                            [r"$\overline{T}_{exo}-25$%", 1.16e-5],
                            [r"$\overline{T}_{exo}+25$%", 3.6e-3],
                            [r"1 pr $\mu m$ H$_2$O", 3.7e-4],
                            [r"100 pr $\mu m$ H$_2$O", 3.6e-4],
                            ["Yung+1988", 0.32], ["Kras. 2000 Model 2", 0.135],
                            ["Kras. 2000 Model 1", 0.016],
                            [r"Kras. 2002 $\odot$ min", 0.055],
                            [r"Kras. 2002 $\odot$ mean", 0.082],
                            [r"Kras. 2002 $\odot$ max", 0.167]],
                           columns=["Name", "Values"])

fig = plt.figure(figsize=(7, 8))
plt.style.use('default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Louis George Caf?'
plt.rcParams['font.monospace'] = 'FreeMono'
plt.rcParams['font.size'] = 19
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22

carray = ["xkcd:faded orange", "xkcd:faded orange", "xkcd:faded orange",
          "xkcd:faded orange", "xkcd:faded orange", "xkcd:faded orange",
          "xkcd:faded orange", "xkcd:faded orange", "xkcd:faded orange",
          "#88bcdf", "#88bcdf", "#88bcdf", "#88bcdf", "#88bcdf", "#88bcdf"]
ax = plt.gca()
ax.set_facecolor("#ededed")
ax.grid(zorder=0, color="white", which="major")
for side in ["top", "bottom", "left", "right"]:
    ax.spines[side].set_visible(False)
ax.barh(list(experiments.index * 2), experiments["Values"], 1.5, color=carray,
        zorder=10)
ax.set_xscale("log")
for i, v in enumerate(experiments["Values"]):
    if v == 1.16e-5:
        m = 5
    else:
        m = -0.1
    ax.text(v+m*v, 2*(i-0.1), str(v), color='black', fontsize=15, zorder=16,
            ha="right")
plt.xlabel("Fractionation Factor")

plt.yticks(range(0, len(experiments)*2, 2), experiments["Name"])
plt.xticks([10**(-5), 10**(-4), 10**(-3), 10**(-2), 10**(-1), 1])

plt.savefig("../Results/VarWaterTemp/Plots/f-plot-pretty-horizontal.png", bbox_inches="tight")