################################################################################
# plot_temp_profs_tradeoffs.jl
# TYPE: Supporting (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: plot the temperature profiles used in the model, one plot each 
# for experiments varying surface, exobase, tropopause temperatures.
# 
# Eryn Cangi
# June 2019
# Currently tested for Julia: 0.7
################################################################################

using PyCall
using PyPlot
using HDF5
using PlotUtils
using LaTeXStrings

function Tpiecewise(z::Float64, Tsurf, Ttropo, Texo, E="")
    #= DO NOT MODIFY! If you want to change the temperature, define a
    new function or select different arguments and pass to Temp(z)

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 "half-Gaussian" function for temperatures 
    altitudes above the tropopause, with a constant lapse rate (1.4K/km) 
    in the lower atmosphere. The tropopause width is allowed to vary
    in certain cases.

    z: altitude above surface in cm
    Tsurf: Surface temperature in K
    Tropo: tropopause tempearture
    Texo: exobase temperature
    E: type of experiment, used for determining if mesopause width will vary 
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of tropopause. It varies unless we're only varying the 
    # exobase temperature.
    if (E=="tropo") || (E=="surf")
        ztropo_bot = (Ttropo-Tsurf)/(lapserate)
        ztropowidth = ztropo - ztropo_bot
    else
        ztropo_bot = (Ttropo-Tsurf)/(lapserate)
        ztropowidth = ztropo - ztropo_bot
    end

    if z >= ztropo  # upper atmosphere
        return Texo - (Texo - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif ztropo > z >= ztropo - ztropowidth  # tropopause
        return Ttropo
    elseif ztropo-ztropowidth > z  # lower atmosphere
        return Tsurf + lapserate*z
    end
end

function get_grad_colors(L, cmap)
    #=
    Generates some colors based on a GRADIENT color map for use in plotting a 
    bunch of lines all at once.
    L: number of colors to generate.
    cmap: color map name

    AVAILABLE MAPS: blues, viridis, pu_or, magma, plasma, inferno

    =#

    colors_dumb = [cgrad(Symbol(cmap))[x] for x in range(0, stop=1, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = red(colors_dumb[i])
        c[i, 2] = green(colors_dumb[i])
        c[i, 3] = blue(colors_dumb[i])
    end
    return c
end

function get_colors(L, cmap)
    #=
    Generates some colors based on a non-gradient color map for use in plotting a 
    bunch of lines all at once.
    L: number of colors to generate.
    cmap: color map name
    =#

    Cmap = get_cmap(Symbol(cmap))
    colors_dumb = [Cmap(x) for x in range(0, stop=1, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = colors_dumb[i][1]
        c[i, 2] = colors_dumb[i][2]
        c[i, 3] = colors_dumb[i][3]
    end
    return c
end

function better_plot_bg(axob)
    axob.set_facecolor("#ededed")
    axob.grid(zorder=0, color="white", which="both")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
end

function make_temp_3panel(A)
    #=
    A: height in km at which to end the altitude grid.
    =#

    Ts = collect(180:10:270)
    Tt = collect(70:10:160)
    Te = collect(125:25:350)
    
    numsurf = length(Ts)
    numtropo = length(Tt)
    numexo = length(Te)

    totallen = length(Ts) + length(Tt) + length(Te)

    alt = collect(0:2e5:A*1e5)

    fig, ax = subplots(1, 3, figsize=(15, 3))
    for axob in ax
        axob.set_xlabel("Temperature (K)")
        axob.set_yticks(collect(0:50:A))
        axob.tick_params(axis="x", labelrotation=45)
        better_plot_bg(axob)
    end
  
    ax[1].set_xticks(collect(110:20:270))
    ax[1].set_xticklabels(collect(110:20:270))
    ax[1].set_ylabel("Altitude (km)")
    ax[1].set_title(L"Varying T$_{surf}$")

    ax[2].set_xticks(collect(70:20:200))
    ax[2].set_xticklabels(collect(70:20:200))
    ax[2].set_title(L"Varying T$_{tropo}$")

    ax[3].set_xticks(collect(125:25:350))
    ax[3].set_xticklabels(collect(125:25:350))
    ax[3].set_title(L"Varying T$_{exo}$")

    for k in range(1, length=totallen)
        csurf = get_colors(numsurf, "plasma")
        if k <= numsurf  # surface
            Tprof = map(h->Tpiecewise(h, Ts[k], 108.0, 205.0, "surf"), alt)
            ax[1].plot(Tprof, alt./1e5, color=csurf[k, :])
        end
        
        csurf = get_colors(numtropo, "plasma")
        if numsurf+1 <= k <= numsurf+numtropo  # mesosphere
            Tprof = map(h->Tpiecewise(h, 216.0, Tt[k-numsurf], 205.0, "tropo"), alt)
            ax[2].plot(Tprof, alt./1e5, color=csurf[k-10, :])
        end
        
        csurf = get_colors(numexo, "plasma")
        if numsurf+numtropo+1 <= k  # exobase
            Tprof = map(h->Tpiecewise(h, 216.0, 108.0, Te[k-(numsurf+numtropo)], "exo"), alt)
            ax[3].plot(Tprof, alt./1e5, color=csurf[k-(numsurf+numtropo), :])
        end
    end
    savefig("../Results/ALL STUDY PLOTS/tradeoff_temp_profiles.png", bbox_inches="tight")
end

rcParams = PyCall.PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 18
rcParams["axes.labelsize"]= 20
rcParams["xtick.labelsize"] = 18
rcParams["ytick.labelsize"] = 18


# =============================================================================

# make_temp_3panel(200)

make_temp_3panel(250)
