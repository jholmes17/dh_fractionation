###################
# plot comparison of water loss results with other studies

# Eryn Cangi
# 23 October 2019
###################

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

studies = ["Yung+ 1988", "Krasnopolsky 2000", "Lammer+ 2003", "Kurokawa+ 2014",
           "Villanueva+ 2015", "Alsaeed+ 2019", "This work"]
dummyind = [0, 1, 2, 3, 4, 5, 6]
water = [[3.4, 3.4], [65, 120], [14, 34], [51, 152], [137, 165], [20, 220],
         [82, 95]]
col = ["xkcd:dodger blue", "xkcd:dodger blue", "xkcd:dodger blue",
       "xkcd:dodger blue", "xkcd:dodger blue", "xkcd:dodger blue", 
       "xkcd:faded orange"]
# col = ["#10007A", "#10007A", "#2F7965", "#2F7965", "#2F7965", "#D59A07",
#        "#10007A"] # "#e16262"
lb = ["1D photochemical modeling", "1D photochemical modeling",
      "Geological observations", "Geological observations",
      "Geological lit + D/H observations", "1D box model",
      "1D photochemical modeling"]

fig, ax = plt.subplots(figsize=(7, 6))
plt.style.use('default')
plt.rc('text', usetex=False)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Louis George Caf?'
plt.rcParams['font.monospace'] = 'FreeMono'
plt.rcParams['font.size'] = 22
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22

# set background/grid
ax.set_facecolor("#ededed")
ax.grid(zorder=0, color="white", which="major")
for side in ["top", "bottom", "left", "right"]:
    ax.spines[side].set_visible(False)

for d, w in zip(dummyind, water):
    ax.plot(w, [d, d], linewidth=15, color=col[d])

p1 = mpatches.Patch(color="#10007A", label="1D photochemical modeling")
p2 = mpatches.Patch(color="#2F7965", label="Observations")
# p3 = mpatches.Patch(color="#e16262", label="Geological lit. + D/H observations")
p3 = mpatches.Patch(color="#D59A07", label="1D box model")

# plt.legend(handles=[p1, p2, p3], fontsize=13, bbox_to_anchor=(0.5, 0.2))

ax.scatter(3.4, 0, marker="|", zorder=10, s=150, color="xkcd:dodger blue")#"#10007A")
plt.yticks(dummyind, labels=studies)
ax.tick_params(axis="both", labelsize=22)
# plt.title("Study comparison of lost water")
plt.xlabel("Water lost (m GEL)", fontsize=22)

plt.savefig("../Results/ALL STUDY PLOTS/waterloss_vs_others.png", bbox_inches="tight")