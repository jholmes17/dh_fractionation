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
# TODO: This has not been updated for the improved temperature function that
# maintains a constant lapse rate and adjusts mesosphere height. 7/2/19
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
    profile_array = [[143, 110, 200], [190, 83, 200], [190, 110, 150], 
    				 [238, 110, 200], [190, 138, 200], [190, 110, 250]]

    titles = [L"$\overline{T}_{surf}-25$%", L"$\overline{T}_{tropo}-25$%", 
              L"$\overline{T}_{exo}-25$%", L"$\overline{T}_{surf}+25$%", 
              L"$\overline{T}_{tropo}+25$%", L"$\overline{T}_{exo}+25$%"]

    rcparams = PyCall.PyDict(matplotlib."rcParams")
    rcparams["font.sans-serif"] = ["Laksaman"]
    rcparams["font.monospace"] = ["FreeMono"]
    rcparams["font.size"] = 22
    rcparams["axes.labelsize"]= 24
    rcparams["xtick.labelsize"] = 22
    rcparams["ytick.labelsize"] = 22

    alt = (0:2e5:200e5)

    fig, ax = subplots(2, 3, sharex=true, sharey=true, figsize=(10,8))
    subplots_adjust(wspace=0.05, hspace=0.4)

    i = 1
    j = 1
    for (el, tit) in zip(profile_array, titles)
        ax[j, i].plot([Tpiecewise(a, el[1], el[2], el[3]) for a in alt], alt/1e5)
        ax[j, i].set_title(tit)
        ax[j, i].set_xticks([100, 150, 200, 250])
        ax[j, i].tick_params(axis="x", labelbottom=false)
        better_plot_bg(ax[j, i])

        if i==1
            ax[j, i].set_ylabel("Altitude (km)")
            ax[j, i].set_yticks([0, 100, 200])
            ax[j, i].set_yticklabels([0,100,200])
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

    ax[1, 1].text(200, 185, L"T_{exo}")
    ax[1, 1].text(115, 75, L"T_{tropo}")
    ax[1, 1].text(155, 10, L"T_{surf}")

    savefig("../Results/VarWaterTemp/Plots/temp_profiles.png")
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

    alt = (0:2e5:200e5)

    fig, ax = subplots(figsize=(4,6))
    better_plot_bg(ax)
    plot([Tpiecewise(a, tprof[1], tprof[2], tprof[3]) for a in alt], alt/1e5)
    title(titletext)

    ax.set_ylabel("Altitude (km)")
    ax.set_yticks([0, 100, 200])
    ax.set_yticklabels([0, 100, 200])
    ax.set_xlabel("Temperature (K)")


    ax.text(tprof[3]-30, 185, L"T_{exo}")
    ax.text(115, 75, L"T_{tropo}")
    ax.text(185, 10, L"T_{surf}")

    savefig(savepath*"/temp_profile_"*titletext*".png", bbox_inches="tight")
end

plot_all_T_profiles()
plot_one_profile([190, 110, 200], "Mean", "../Results/VarWaterTemp/Plots")