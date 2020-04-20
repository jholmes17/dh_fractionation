################################################################################
# plot_temp_profs_tradeoffs.jl
# TYPE: Supporting (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: plot the temperature profiles used in the model, one plot each 
# for experiments varying surface, exobase, tropopause temperatures.
# 
# Eryn Cangi
# June 2019
# Last edited: 3 Janauary 2020
# Testing status: YELLOW (changes made, but untested)
# Currently tested for Julia: 0.7
################################################################################

using PyCall
using PyPlot
using HDF5
using PlotUtils
using LaTeXStrings
using Analysis

include("PARAMETERS.jl")


function make_temp_3panel(base)
    #=
    base: main results folder (general)
    =#

    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    Ts = collect(150:10:270)
    Tt = collect(lowTt:10:hiTt)
    Te = collect(150:25:350)
    
    numsurf = length(Ts)
    numtropo = length(Tt)
    numexo = length(Te)
    totallen = length(Ts) + length(Tt) + length(Te)

    fig, ax = subplots(1, 3, figsize=(15, 3))
    for axob in ax
        axob.set_xlabel("Temperature (K)")
        axob.set_yticks(collect(0:50:zmax/1e5))
        axob.tick_params(axis="x", labelrotation=45)
        plot_bg(axob)
    end
  
    # Surface axis
    ax[1].set_xticks(collect(meanTt:20:hiTs))
    ax[1].set_xticklabels(collect(meanTt:20:hiTs))
    ax[1].set_ylabel("Altitude (km)")
    ax[1].set_title(L"Varying T$_{surf}$")

    # Tropopause axis
    ax[2].set_xticks(collect(lowTt:20:220))
    ax[2].set_xticklabels(collect(lowTt:20:220))
    ax[2].set_title(L"Varying T$_{tropo}$")

    # Exobase axis 
    ax[3].set_xticks(collect(125:25:350))
    ax[3].set_xticklabels(collect(125:25:350))
    ax[3].set_title(L"Varying T$_{exo}$")

    for k in range(1, length=totallen)
        cols = get_colors(numsurf, "plasma")
        if k <= numsurf  # surface
            Tprof = map(h->Tpiecewise(h, Ts[k], meanTt, meanTe, "surf"), alt)
            ax[1].plot(Tprof, alt./1e5, color=cols[k, :])
        end
        
        cols = get_colors(numtropo, "plasma")
        if numsurf+1 <= k <= numsurf+numtropo  # mesosphere
            Tprof = map(h->Tpiecewise(h, meanTs, Tt[k-numsurf], meanTe, "tropo"), alt)
            ax[2].plot(Tprof, alt./1e5, color=cols[k-numsurf, :])
        end
        
        cols = get_colors(numexo, "plasma")
        if numsurf+numtropo+1 <= k  # exobase
            Tprof = map(h->Tpiecewise(h, meanTs, meanTt, Te[k-(numsurf+numtropo)], "exo"), alt)
            ax[3].plot(Tprof, alt./1e5, color=cols[k-(numsurf+numtropo), :])
        end
    end
    savefig(base*"ALL STUDY PLOTS/tradeoff_temp_profiles.png", bbox_inches="tight")
    savefig(base*"TradeoffPlots/Main Results/tradeoff_temp_profiles.png", bbox_inches="tight")
end


# =============================================================================

# make_temp_3panel(200)
base = "/home/emc/GDrive-CU/Research/Results/"
make_temp_3panel(base)
