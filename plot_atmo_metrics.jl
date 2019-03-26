################################################################################
# Plot metrics of atmosphere like CO/O2, H2, O3.
# currently assumes a bunch of sub-experiment folders.
#
# Eryn Cangi 
# 11 March 2019
# Currently tested for Julia: 0.7
################################################################################

# modules ---------- #
using PyPlot
using HDF5
using LaTeXStrings
# ------------------ #

# utility functions ------------------------------------------------------------
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
#-------------------------------------------------------------------------------

numexp = 9  # number of experiments for which to plot values
# base = "/data/GDrive-CU/Research/Results/VarWaterTemp/"
base = "/home/emc/GDrive-CU/Research/Results/VarWaterTemp/"
expfolders = ["temp_192_83_199", "temp_192_110_149", "temp_192_110_199", 
              "temp_240_110_199", "temp_192_138_199", "temp_192_110_249", 
              "water_1e-3", "water_1e-4", "water_1e-5"]
explbls = ["low tropo", "low exo", "mean", "high surface", "high tropo", 
           "high exo", "water 1e-3", "water 1e-4", "water 1e-5"]

# to plot Mike's results
# numexp = 1
# base = "/data/GDrive-CU/Research/"
# expfolders = ["chaffin_natgeo_mars_photochemistry"]
# explbls = ["H model only"]

# CO/O2 ratio ------------------------------------------------------------------
oxy_vals = Array{Float64}(numexp,99)
oxyvals_surf = Array{Float64}(numexp)
H2ppm = Array{Float64}(numexp)
O3ppm = Array{Float64}(numexp)

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


# Figures ----------------------------------------------------------------------
fig = figure(figsize=(6,4))
ax = gca()
bar(collect(1:numexp), oxyvals_surf, zorder=2)
ax[:axhline](0.6, color="black", label="Observed [Nair 1994")
grid(color="gainsboro", zorder=1)
xlabel("Experimental runs of model")
ylabel(L"CO/O_2")
ax[:set_xticks](1:numexp)
ax[:set_xticklabels](explbls[1:numexp], rotation=30)
title("Oxidation in equilibrium")
legend()
tight_layout()
savefig(base*"CO-O2-ratio.png")
show()

fig2 = figure(figsize=(6,4))
ax = gca()
bar(collect(1:numexp), H2ppm, zorder=2)
grid(color="gainsboro", zorder=1)
ax[:axhline](15, label="Krasnopolsky 2001", color="black")
xlabel("Experimental runs of model")
ylabel(L"H_2"*" (ppm)")
ax[:set_xticks](1:numexp)
ax[:set_xticklabels](explbls[1:numexp], rotation=30)
title(L"H_2"*" concentration in equilibrium")
legend()
tight_layout()
savefig(base*"H2-MR.png")
show()

fig3 = figure(figsize=(6,4))
ax = gca()
bar(collect(1:numexp), O3ppm, zorder=2)
grid(color="gainsboro", zorder=1)
ax[:axhline](10^10, label="Lefevre 2004", color="black")
xlabel("Experimental runs of model")
ylabel(L"O_3"*" (column abundance)")
ax[:set_xticks](1:numexp)
ax[:set_xticklabels](explbls[1:numexp], rotation=30)
yscale("log")
ylim(10^9, 3*10^10)
title(L"O_3"*" concentration in equilibrium")
legend()
tight_layout()
savefig(base*"O3.png")
show()