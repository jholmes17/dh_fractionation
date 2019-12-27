################################################################################
# plot_atmo_metrics.jl
# TYPE: MAIN (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: Make a 3-panel verticla plot of atmospheric metrics CO/O2, H2, O3.
# currently assumes a bunch of sub-experiment folders. Used to generate this plot
# for poster in 2019 for various conferences including 1st MAVEN PSG of the year.
#
# Eryn Cangi 
# 11 March 2019
# Currently tested for Julia: 0.7
################################################################################

# modules ======================================================================
using PyPlot
using HDF5
using LaTeXStrings
using PyCall

# utility functions ============================================================
function get_ncurrent(readfile)
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],
                                                          length(alt)-2)
    end
    n_current
end

const alt = Array((0:2e5:200e5))

numexp = 8  # number of experiments for which to plot values
# base = "/data/GDrive-CU/Research/Results/VarWaterTemp/"
base = "/home/emc/GDrive-CU/Research/Results/VarWaterTemp/"
expfolders = ["temp_192_110_199", "temp_192_83_199", "temp_192_110_149", 
              "temp_240_110_199", "temp_192_138_199", "temp_192_110_249", 
              "water_1e-3", "water_1e-5"]
explbls = ["Mean", "Low tropo", "Low exo", "High surface", "High tropo", 
           "High exo", "High water", "Low water"]

# to plot Mike's results
# numexp = 1
# base = "/data/GDrive-CU/Research/"
# expfolders = ["chaffin_natgeo_mars_photochemistry"]
# explbls = ["H model only"]

# Get data =====================================================================
oxy_vals = Array{Float64}(undef, numexp,99)
oxyvals_surf = Array{Float64}(undef, numexp)
H2ppm = Array{Float64}(undef, numexp)
O3ppm = Array{Float64}(undef, numexp)

for i in 1:1:numexp
    f = expfolders[i]
    # in following line get rid of _D_"*f*" if doing Mike's results
    alldata = h5read(base*f*"/converged_standardwater_D_"*f*".h5", 
                     "n_current/n_current_mat")
    ntot = sum(alldata)
    atmo = get_ncurrent(base*f*"/converged_standardwater_D_"*f*".h5")

    # CO/O2 ratio --------------------------------------------------------------
    COdata = atmo[:CO]
    O2data = atmo[:O2]
    oxy_vals[i, :] = COdata./O2data;
    oxyvals_surf[i] = signif(oxy_vals[i, :][1], 2)
    println("CO/O2, "*explbls[i]*": ", signif(oxy_vals[i, :][1], 2))
    println()

    # H2 -----------------------------------------------------------------------
    H2data = atmo[:H2]
    H2ppm[i] = signif((sum(H2data)/ntot)/1e-6, 3)
    println("H2 ppm, "*explbls[i]*": ", signif((sum(H2data)/ntot)/1e-6, 3))
    println()

    # O3 -----------------------------------------------------------------------
    O3data = atmo[:O3]
    O3ppm[i] = signif(sum(O3data), 3)
    println("O3 column, "*explbls[i]*": ", signif(sum(O3data), 3))
    println()
end


# Figures ======================================================================
# make plots pretty
rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 22
rcParams["axes.labelsize"]= 24
rcParams["xtick.labelsize"] = 22
rcParams["ytick.labelsize"] = 22


# Make a stacked version of all interesting metrics
fig, axs = subplots(3, 1, sharex=true, figsize=(11,21))
subplots_adjust(hspace=0, bottom=0.15)
# CO/O2
axs[1][:bar](collect(1:numexp), oxyvals_surf, zorder=2, 
             color="xkcd:faded orange")
axs[1][:axhline](0.6, color="black", label="Owen 1977")
axs[1][:grid](color="gainsboro", zorder=1)
axs[1][:set_ylabel](L"CO/O_2")
axs[1][:legend](loc="upper left", bbox_to_anchor=(0.55, 0.99))
# H2
axs[2][:bar](collect(1:numexp), H2ppm, zorder=2, color="xkcd:faded orange")
axs[2][:grid](color="gainsboro", zorder=1)
axs[2][:axhline](15, label="Krasnopolsky & Feldman 2001", color="black")
axs[2][:fill_between](collect(0:numexp+1), 10, 20, alpha=0.3, color="black", 
                      zorder=3)
axs[2][:set_ylabel](L"H_2"*" (ppm)")
axs[2][:legend](fontsize=21)
# O3 Ozone
axs[3][:bar](collect(1:numexp), O3ppm, zorder=2, color="xkcd:faded orange")
axs[3][:grid](color="gainsboro", zorder=1)
axs[3][:axhline](0.75*10^10, label="Lefèvre 2004", color="black")
axs[3][:fill_between](collect(0:numexp+1), 0.5*10^10, 10^10, alpha=0.3, 
                      color="black", zorder=3)
axs[3][:set_ylabel](L"O_3"*" (cm"*L"^{-3}"*")")
axs[3][:set_yscale]("log")
axs[3][:set_xlim](0.5,8.5)
axs[3][:set_ylim](10^9, 3*10^10)
axs[3][:set_xticks](1:numexp)
axs[3][:set_xticklabels](explbls[1:numexp], rotation=30, ha="right")
axs[3][:legend]()
savefig(base*"metrics_stacked.png")
show()

# Other: individual plots ======================================================
# vertical profile oxidation
# fig, ax = subplots(1, numexp, figsize=(30,6))
# # do not change the following and do not use tight_layout(), it conflicts
# subplots_adjust(wspace=0, top=0.9, bottom=0.15)
# #grid(color="gainsboro", zorder=1)
# suptitle(L"CO/O_2")#, y=1.05)

# for i in 1:8
#     ax[i][:plot](oxy_vals[i, :], (2:2:198), zorder=2, color="xkcd:faded orange")
#     ax[i][:axvline](0.6, color="#88bcdf", label="Observed=0.6\n[Nair 1994]")
#     ax[i][:set_xlabel](explbls[i])
#     ax[i][:set_xscale]("log")
# end
# ax[1][:set_ylabel]("Altitude (km)")

# savefig(base*"CO-O2-ratio-vert.png")
# show()

# surface oxidation
# fig = figure(figsize=(11,7))
# ax = gca()
# bar(collect(1:numexp), oxyvals_surf, zorder=2, color="xkcd:faded orange")
# ax[:axhline](0.6, color="black", label="Observed [Nair 1994]")
# grid(color="gainsboro", zorder=1)
# ylabel(L"CO/O_2")
# ax[:set_xticks](1:numexp)
# ax[:set_xticklabels](explbls[1:numexp], rotation=30, ha="right")
# #title("Oxidation in equilibrium", fontsize=36)
# legend(loc="upper left", bbox_to_anchor=(0.37, 0.95))
# tight_layout()
# savefig(base*"CO-O2-ratio.png")
# show()

# H2 population
# fig2 = figure(figsize=(11,7))
# ax = gca()
# bar(collect(1:numexp), H2ppm, zorder=2, color="xkcd:faded orange")
# grid(color="gainsboro", zorder=1)
# ax[:axhline](15, label="Krasnopolsky 2001", color="black")
# fill_between(collect(0:numexp+1), 10, 20, alpha=0.3, color="black", zorder=3)
# ylabel(L"H_2"*" (ppm)")
# ax[:set_xticks](1:numexp)
# ax[:set_xticklabels](explbls[1:numexp], rotation=30, ha="right")
# #title(L"H_2"*" concentration in equilibrium", fontsize=36)
# legend()
# xlim(0, 9)
# tight_layout()
# savefig(base*"H2-MR.png")
# show()

# Ozone. Note: the Lefevre data is from his 3D modeling paper, Figure 2d. 
# fig3 = figure(figsize=(11,7))
# ax = gca()
# bar(collect(1:numexp), O3ppm, zorder=2, color="xkcd:faded orange")
# grid(color="gainsboro", zorder=1)
# ax[:axhline](0.75*10^10, label="Lefèvre 2004", color="black")
# fill_between(collect(0:numexp+1), 0.5*10^10, 10^10, alpha=0.3, color="black", zorder=3)
# ylabel(L"O_3"*" (column abundance)")
# ax[:set_xticks](1:numexp)
# ax[:set_xticklabels](explbls[1:numexp], rotation=30, ha="right")
# yscale("log")
# xlim(0,9)
# ylim(10^9, 3*10^10)
# #title(L"O_3"*" concentration in equilibrium", fontsize=36)
# legend()
# tight_layout()
# savefig(base*"O3.png")
# show()
