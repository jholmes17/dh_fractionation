################################################################################
# plot_temp_profs.jl
# TYPE: Supporting (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: plot the temperature profiles used in the model, either
# all in a 6-panel plot or individual profiles.
# 
# Eryn Cangi
# May 2019
# Currently tested for Julia: 0.7
################################################################################

using PyCall
using PyPlot
using HDF5
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

function better_plot_bg(axob)
    axob.set_facecolor("#ededed")
    axob.grid(zorder=0, color="white", which="both")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
end

function plot_all_T_profiles()
    #=
    Plots a 6-panel figure of temperature profiles, the ones where there is 
    some perturbance from mean.
    =#
    profile_array = [[162, 108, 205], [216, 81, 205], [216, 108, 158], 
    				 [270, 108, 205], [216, 135, 205], [216, 108, 264]]

    titles = [L"$T_{surf}=0.75\overline{T}_{surf}$", L"$T_{tropo}=0.75\overline{T}_{tropo}$", 
              L"$T_{exo}=0.75\overline{T}_{exo}$", L"$T_{surf}=1.25\overline{T}_{surf}$", 
              L"$T_{tropo}=1.25\overline{T}_{tropo}$", L"$T_{exo}=1.25\overline{T}_{exo}$"]

    rcparams = PyCall.PyDict(matplotlib."rcParams")
    rcparams["font.sans-serif"] = ["Louis George Caf?"]
    rcparams["font.monospace"] = ["FreeMono"]
    rcparams["font.size"] = 22
    rcparams["axes.labelsize"]= 24
    rcparams["xtick.labelsize"] = 22
    rcparams["ytick.labelsize"] = 22

    alt = (0:2e5:250e5)

    fig, ax = subplots(2, 3, sharex=true, sharey=true, figsize=(10,8))
    subplots_adjust(wspace=0.05, hspace=0.3)

    i = 1
    j = 1
    modcol = "navy"
    meancol = "xkcd:cerulean blue"
    for (el, tit) in zip(profile_array, titles)
        ax[j, i].plot([Tpiecewise(a, 216.0, 108.0, 205.0) for a in alt], 
                      alt/1e5, color=meancol, linestyle="--", label="Global mean")
        ax[j, i].plot([Tpiecewise(a, el[1], el[2], el[3]) for a in alt], 
                       alt/1e5, color=modcol, label="Modified") 
        ax[j, i].set_title(tit, fontsize=20)
        ax[j, i].set_xticks([100, 150, 200, 250])
        ax[j, i].tick_params(axis="x", labelbottom=false)
        better_plot_bg(ax[j, i])

        if i==1
            ax[j, i].set_ylabel("Altitude (km)")
            ax[j, i].set_yticks([0, 50, 100, 150, 200, 250])
            ax[j, i].set_yticklabels([0,50,100,150,200,250])
            if j==2
                ax[j, i].set_xlabel("Temperature (K)")
                ax[j, i].tick_params(axis="x", labelbottom=true)
            end
        end
        
        i+=1
        if i == 4
            i = 1
            j = 2
        else
            continue
        end
    end

    ax[1, 1].text(213, 220, L"T_{exo}")
    ax[1, 1].text(115, 80, L"T_{tropo}")
    ax[1, 1].text(80, -2, L"T_{surf-}")
    ax[1, 1].plot(205, 250, marker="o", color=modcol)
    ax[1, 1].plot(108, 80, marker="o", color=modcol)
    ax[1, 1].plot(162, 0, marker="o", color=modcol)
    legend(bbox_to_anchor=(1.05, -0.05))

    savefig("../Results/ALL STUDY PLOTS/temp_profiles.png", bbox_inches="tight")
end

function plot_one_profile(tprof, titletext, savepath)
    #=
    tprof: a list of the temps [T_surface, T_tropo, T_exo]
    titletext: exactly what it sounds like
    savepath: where to save the figure

    This function also appears within the code of converge_new_file, but this 
    function can be called if for some reason a profile wasn't generated.
    =#
    # make plots pretty
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    alt = (0:2e5:250e5)

    fig, ax = subplots(figsize=(4,6))
    better_plot_bg(ax)
    plot([Tpiecewise(a, tprof[1], tprof[2], tprof[3]) for a in alt], alt/1e5, color="xkcd:cerulean blue")
    title(titletext)

    ax.set_ylabel("Altitude (km)")
    ax.set_xticks([100,125,150,175,200,225  ])
    ax.set_yticks([0, 50, 100, 150, 200, 250])
    ax.set_yticklabels([0,50,100,150,200,250])
    ax.set_xlabel("Temperature (K)")

    ax.text(tprof[3]-25, 230, L"\overline{T}_{exo}")
    ax.text(110, 85, L"\overline{T}_{tropo}")
    ax.text(175, 0, L"\overline{T}_{surf}")

    savefig(savepath*"/temp_profile_"*titletext*".png", bbox_inches="tight")
end

# plot_all_T_profiles()
plot_one_profile([216, 108, 205], "Global mean temperature", "../Results/ALL STUDY PLOTS/")