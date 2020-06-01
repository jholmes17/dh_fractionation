################################################################################
# plot_f_results.jl
# TYPE: MAIN (analysis)
# WHICH: Equilibrium experiments
# DESCRIPTION: easily calculate f for a given converged file. 
#
# Eryn Cangi
# 31 May 2019
# Last edited: 21 April 2020
# Currently tested for Julia: 1.4.1
################################################################################

using PyPlot
using PyCall
using HDF5
using Printf
using LaTeXStrings
using Analysis
using DataFrames

include("PARAMETERS.jl")

global rcParams = PyCall.PyDict(matplotlib."rcParams")

# Functions to calculate f for a converged Yung case experiment ================
function yungTemp(z::Float64)
    # Temperature profile from Kong & McElroy 1977
    kongtemps = Array{Float64}([220., 214., 208., 202.2, 196.6, 191., 185.6, 180.2, 174.8, 169.4, 164., 161.6,
                                159.2, 156.8, 154.4, 152., 150.6, 149.2, 147.8, 146.4, 145., 144., 143., 142.,
                                141., 140., 139.4, 138.8, 138.2, 137.6, 137., 137., 137., 137., 137., 137.,
                                137., 137., 137., 137., 137., 137.6, 138.2, 138.8, 139.4, 140., 145., 150.,
                                155., 160., 165., 171.4, 177.8, 184.2, 190.6, 197., 203.4, 209.8, 216.2, 222.6,
                                229., 235.6, 242.2, 248.8, 255.4, 262., 268.8, 275.6, 282.4, 289.2, 296., 301.8,
                                307.6, 313.4, 319.2, 325., 328.2, 331.4, 334.6, 337.8, 341., 343.2, 345.4, 347.6,
                                349.8, 352., 353.4, 354.8, 356.2, 357.6, 359., 359.6, 360.2, 360.8, 361.4, 362.,
                                362.4, 362.8, 363.2, 363.6, 364., 364.2, 364.4, 364.6, 364.8, 365., 365., 365.,
                                365., 365., 365., 365., 365., 365., 365., 365., 365., 365., 365., 365.,
                                365., 365., 365., 365., 365., 365., ])
    i = Int(z / 2e5) + 1
    return kongtemps[i]
end

function speciesbcs_YUNG(species, oflux)
    #=
    Defines boundary conditions for a variable O flux and temperature set.
    species: symbol; the species to get the boundary condition.
    oflux: a number such as 1.2e8
    temps: an array of surface temp, tropopause temp, exobase temp. 
    =#
    H2Osat = map(x->Psat(x), map(yungTemp, alt))
    HDOsat = map(x->Psat_HDO(x), map(yungTemp, alt))
    H_effusion_velocity = effusion_velocity(yungTemp(zmax),1.0)
    H2_effusion_velocity = effusion_velocity(yungTemp(zmax),2.0)
    D_effusion_velocity = effusion_velocity(yungTemp(zmax),2.0)
    HD_effusion_velocity = effusion_velocity(yungTemp(zmax),3.0)

    speciesbclist=Dict(
                    :CO2=>["n" 2.1e17; "f" 0.],
                    :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                    :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                    :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                    :HDO=>["n" HDOsat[1]; "f" 0.],
                    :O=>["f" 0.; "f" oflux],
                    :H2=>["f" 0.; "v" H2_effusion_velocity],
                    :HD=>["f" 0.; "v" HD_effusion_velocity],
                    :H=>["f" 0.; "v" H_effusion_velocity],
                    :D=>["f" 0.; "v" D_effusion_velocity],
                   );
    get(speciesbclist,
        species,
        ["f" 0.; "f" 0.])
end

function get_H_flux_yung(readfile, oflux)
    #=
    Produces an array of H flux out of the top of the atmosphere, for a given
    O flux and temperature array.
    readfile: a converged simulation, type .h5
    oflux: a value such as 1.2e8 representing the O flux boundary condition
    temps: an array of surface temp, tropopause temp, exobase temp.  
    =#

    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (H population @ 200 km) * (D flux rate)
    # and (H2 population @ 200 km) * 2*(H2 flux rate) and 
    # (HD population @ 200 km) * (HD flux rate)
    Hfluxes = (n_current[:H][end]*speciesbcs_YUNG(:H, oflux)[2,2]
                  + 2*n_current[:H2][end]*speciesbcs_YUNG(:H2, oflux)[2,2]
                  + n_current[:HD][end]*speciesbcs_YUNG(:HD, oflux)[2,2])
    return Hfluxes
end

function get_D_flux_yung(readfile, oflux)
    #=
    Gets value of D flux out of the atmosphere in a converged file.
    =#
    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (D population @ 200 km) * (D flux rate)
    # and (HD population @ 200 km) * (HD flux rate)
    Dfluxes = (n_current[:D][end]*speciesbcs_YUNG(:D, oflux)[2,2]
                  + n_current[:HD][end]*speciesbcs_YUNG(:HD, oflux)[2,2])
    return Dfluxes
end

function calculate_f_yung(thefile, Oflux)

    ncur = get_ncurrent(thefile)
    Hf = get_H_flux_yung(thefile, Oflux)
    Df = get_D_flux_yung(thefile, Oflux)
    return 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1])
end

function do_it_for_yung_case()
    n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
    calculate_f_yung(results_dir*"Yung-With-Old-Water/Case 1/2019 redo for assurance/converged_yung_case1.h5", 1.2e8)
end


# Main =========================================================================
function plot_results_barchart(escape_type)
    #=
    escape_type: "thermal" or "both". Makes a plot of fractionation factor by 
                 experiment for escape_tpe of atmospheric escape (thermal, or 
                 thermal + nonthermal).
    =#
    n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
    Oflux = 1.2e8# whatever Oflux is, in format "1.2e8" cm^-2 s^-1
    base = results_dir*"VarWaterTemp/"

    # meanTs = 216.0
    # downTs = 162.0
    # upTs = 270.0
    # meanTt = 108.0
    # downTt = 81.0
    # upTt = 135.0
    # meanTe = 205.0
    # downTe = 158.0
    # upTe = 264.0

    flist_mine = [["temp_$(Int(meanTs))_$(Int(meanTt))_$(Int(meanTe))", "Mean profile", 0],
                 ["temp_$(Int(downTs))_$(Int(meanTt))_$(Int(meanTe))", L"$0.75\overline{T}_{surf}$", 0],
                 ["temp_$(Int(upTs))_$(Int(meanTt))_$(Int(meanTe))", L"$1.25\overline{T}_{surf}$", 0],
                 ["temp_$(Int(meanTs))_$(Int(downTt))_$(Int(meanTe))", L"$0.75\overline{T}_{tropo}$", 0],
                 ["temp_$(Int(meanTs))_$(Int(upTt))_$(Int(meanTe))", L"$1.25\overline{T}_{tropo}$", 0],
                 ["temp_$(Int(meanTs))_$(Int(meanTt))_$(Int(downTe))", L"$0.75\overline{T}_{exo}$", 0],
                 ["temp_$(Int(meanTs))_$(Int(meanTt))_$(Int(upTe))", L"$1.25\overline{T}_{exo}$", 0],
                 ["water_2.1e-5", L"1 pr $\mu m$ H$_2$O", 0],
                 # ["water_2.1e-5", L"1 pr $\mu m$ H$_2$O", 0],
                 ["water_2.18e-3", L"25 pr $\mu m$ H$_2$O", 0],
                 ["past", "Yung+1988", 0.32],
                 ["past", L"Kras. 2002 $\odot$ min", 0.055],
                 ["past", L"Kras. 2002 $\odot$ mean", 0.082],
                 ["past", L"Kras. 2002 $\odot$ max", 0.167], 
                 ["past", "Kras. 2000 Model 2", 0.135],
                 ["past", "Kras. 2000 Model 1", 0.016]
                 ]

    # the next line is witchcraft that makes the above array make sense. stupid julia.
    # i hate the way Julia handles array shapes. 
    flist_mine = permutedims(reshape(hcat(flist_mine...), (length(flist_mine[1]), length(flist_mine))))
    
    for f in readdir(base)
        elems = split(f, "_")
        # FNext = f
        if elems[1] == "temp"
            t = [parse(Float64, elems[2]), parse(Float64, elems[3]), parse(Float64, elems[4])]
        elseif elems[1]=="water"
            t = [meanTs, meanTt, meanTe]
        else
            continue
        end
        filetouse = base*f*"/converged_"*f*".h5"

        for row in 1:1:size(flist_mine)[1]
            if flist_mine[row, 1] == f  # ignore the past studies
                flist_mine[row, 3] = calculate_f(filetouse, escape_type, t, Oflux)
            end
        end
    end

    standout = "xkcd:faded orange"
    oldstuff = "#bdd2db" #  "#88bcdf"
    carray_mine = [standout, standout, standout, standout, standout, standout, 
                   standout, standout, standout, oldstuff, oldstuff, oldstuff, 
                   oldstuff, oldstuff, oldstuff]

    fig2, ax = subplots(figsize=(7,7))
    better_plots(ax)
    ax.barh(collect(1:1:size(flist_mine)[1]), flist_mine[:, 3], 0.9, color=carray_mine, zorder=10)
    ax.set_xlabel("Fractionation Factor", fontsize=20)
    ax.set_xscale("log")
    ax.set_yticks(collect(1:1:size(flist_mine)[1]))
    ax.set_yticklabels(flist_mine[:, 2])
    ax.tick_params(which="both", labelsize=18)  # for some reason this has to be here because editing rcParams in a function doesn't work if you do it more than once.
    for (i, v) in enumerate(flist_mine[:, 3])
        if v <= 2e-5
            m = 4
        else
            m = -0.05
        end

        if escape_type=="thermal"
            if flist_mine[i, 1] != "past"
                printvalstyle = " $(@sprintf("%.1e", v))"
            else
                printvalstyle = "$(round(v, digits=3))"
            end
        elseif escape_type=="both"
            printvalstyle = "$(round(v, digits=3))"
        else
            println("Invalid escape type: $(escape_type)")
        end
        ax.set_xticks([10^(-5), 10^(-4), 10^(-3), 10^(-2), 10^(-1), 1])
        ax.text(v+m*v, i-0.1, printvalstyle, color="black", fontsize=15,zorder=16,
                ha="right")
        if escape_type=="thermal"
            plt.title("Thermal escape only")
        else
            plt.title("Thermal + nonthermal escape")
        end
    end

    savefig(results_dir*"ALL STUDY PLOTS/f-results-plot-$(escape_type).png", bbox_inches="tight")
end

function plot_results_caltech_together(base)
    #=
    Plots the thermal and thermal+nonthermal results together in the same plot.

    base: the folder of results
    =#
    n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
    Oflux = 1.2e8# whatever Oflux is, in format "1.2e8" cm^-2 s^-1
    

    # get high and low f in surfaces
    f_surf = Array{Float64}(undef, 3, 4) #array to store f
    f_surf[:, 1] = [lowTs, meanTs, hiTs]
    i = 1
    for Ts in f_surf[:, 1]
        # construct file name
        F = "temp_$(Int(Ts))_$(Int(meanTt))_$(Int(meanTe))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_surf[i, 2] = calculate_f(filetouse, "thermal", [Ts, meanTt, meanTe], Oflux)
        f_surf[i, 3] = calculate_f(filetouse, "both", [Ts, meanTt, meanTe], Oflux)
        f_surf[i, 4] = calculate_f(filetouse, "nonthermal", [Ts, meanTt, meanTe], Oflux)
        i += 1
    end

    # get high and low f in tropopause
    f_tropo = Array{Float64}(undef, 3, 4) #array to store f
    f_tropo[:, 1] = [lowTt, meanTt, hiTt]
    i = 1
    for Tt in f_tropo[:, 1]
        # construct file name
        F = "temp_$(Int(meanTs))_$(Int(Tt))_$(Int(meanTe))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_tropo[i, 2] = calculate_f(filetouse, "thermal", [meanTs, Tt, meanTe], Oflux)
        f_tropo[i, 3] = calculate_f(filetouse, "both", [meanTs, Tt, meanTe], Oflux)
        f_tropo[i, 4] = calculate_f(filetouse, "nonthermal", [meanTs, Tt, meanTe], Oflux)
        i += 1
    end

    # get high and low f in exobase
    f_exo = Array{Float64}(undef, 3, 4) #array to store f
    f_exo[:, 1] = [lowTe, meanTe, hiTe]
    i = 1
    for Te in f_exo[:, 1]
        # construct file name
        F = "temp_$(Int(meanTs))_$(Int(meanTt))_$(Int(Te))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_exo[i, 2] = calculate_f(filetouse, "thermal", [meanTs, meanTt, Te], Oflux)
        f_exo[i, 3] = calculate_f(filetouse, "both", [meanTs, meanTt, Te], Oflux)
        f_exo[i, 4] = calculate_f(filetouse, "nonthermal", [meanTs, meanTt, Te], Oflux)
        i += 1
    end

    # get high and low f in water
    f_water = Array{Float64}(undef, 5, 4) #array to store f
    waterfolders = search_subfolders(base, r"water_\d.+")#r"water_[0-9]+\.[0-9]+e-[0-9]")
    i = 1
    for w in waterfolders

        # get the experiment name
        waterexp = match(r"water_\d.+", w).match
        watermr = match(r"\d.+", waterexp).match
        # construct file name
        filetouse = w * "/converged_"*waterexp*".h5"

        # calculate f
        f_water[i, 1] = parse(Float64, watermr)
        f_water[i, 2] = calculate_f(filetouse, "thermal", [meanTs, meanTt, meanTe], Oflux)
        f_water[i, 3] = calculate_f(filetouse, "both", [meanTs, meanTt, meanTe], Oflux)
        f_water[i, 4] = calculate_f(filetouse, "nonthermal", [meanTs, meanTt, meanTe], Oflux)
        i += 1
    end

    skey = L"T$_{surface}$"
    tkey = L"T$_{tropopause}$"
    ekey = L"T$_{exobase}$"

    # iterate through all the detailed files and make an array of fractionation factors by experiment
    # the mean case file is needed
    mn = base * "temp_$(meanTsint)_$(meanTtint)_$(meanTeint)/converged_temp_$(meanTsint)_$(meanTtint)_$(meanTeint).h5"

    # the f calculations for thermal 
    f_mean_thermal = calculate_f(mn, "thermal", [meanTs, meanTt, meanTe], Oflux)
    f_mean_both = calculate_f(mn, "both", [meanTs, meanTt, meanTe], Oflux)
    f_mean_nonthermal = calculate_f(mn, "nonthermal", [meanTs, meanTt, meanTe], Oflux)

    f_thermal = DataFrame(Exp=["Surface", "Tropopause", "Exobase", "Water", ""],
                          Min=[minimum(f_surf[:, 2]), minimum(f_tropo[:, 2]), 
                               minimum(f_exo[:, 2]), minimum(f_water[:, 2]), 1],
                          Max=[maximum(f_surf[:, 2]), maximum(f_tropo[:, 2]), 
                               maximum(f_exo[:, 2]), maximum(f_water[:, 2]), 1])

    f_both = DataFrame(Exp=["Surface", "Tropopause", "Exobase", "Water", ""],
                       Min=[minimum(f_surf[:, 3]), minimum(f_tropo[:, 3]), 
                             minimum(f_exo[:, 3]), minimum(f_water[:, 3]), 1],
                       Max=[maximum(f_surf[:, 3]), maximum(f_tropo[:, 3]), 
                            maximum(f_exo[:, 3]), maximum(f_water[:, 3]), 1])

    f_nonthermal = DataFrame(Exp=["Surface", "Tropopause", "Exobase", "Water", ""],
                             Min=[minimum(f_surf[:, 4]), minimum(f_tropo[:, 4]), 
                                  minimum(f_exo[:, 4]), minimum(f_water[:, 4]), 1],
                             Max=[maximum(f_surf[:, 4]), maximum(f_tropo[:, 4]), 
                                  maximum(f_exo[:, 4]), maximum(f_water[:, 4]), 1])
 
    # Print the results for use in a table in the paper    
    println("thermal")
    println("f mean: $(f_mean_thermal)")
    println(f_thermal)
    println()
    println("thermal + nonthermal")
    println("f mean: $(f_mean_both)")
    println(f_both)
    println()
    println("nonthermal")
    println("f mean: $(f_mean_nonthermal)")
    println(f_nonthermal)

    toplot_others = Dict("Clarke+ 2019"=>[0.016, 0.047],
                   "Krasnopolsky 2002"=>[0.055, 0.082, 0.167],
                   "Krasnopolsky 2000"=>[0.016, 0.135],
                   "Yung+ 1988"=>[0.32]
                  )


    # PLOT =====================================================================
    fig, ax = subplots(figsize=(12,8))

    # STYLE SETTINGS -----------------------------------------------------------
    for side in ["left", "right"]
        ax.spines[side].set_visible(false)
    end
    ax.tick_params(which="both", labelleft=false, left=false, labeltop=true, top=true)

    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    # Which things to plot 
    plot_order_mine = ["Surface", "Tropopause", "Exobase", "Water", ""]
    plot_order_others = ["Clarke+ 2019", "Krasnopolsky 2002", 
                         "Krasnopolsky 2000", "Yung+ 1988"]
    dummyind_mine = [1, 2, 3, 4, 5]
    dummyind_others = [6, 7, 8, 9]

    myorange = "xkcd:faded orange"
    colors = ["#A35B24", myorange, "#FFAB6B", "#1493A3",
            #myorange, myorange, myorange, myorange, 
             "white", "#3c6843", "#682f51", "#99801d", "#485468"]
    scatcolors = ["#6B3C18", "#B8713B", "#E79B61", "#0e707c", "#fff4e0", 
                  "#3c6843", "#682f51", "#99801d", "#485468"]

    m = "D"
    lw = 5
    msize = 60 
    nameha = "right"
    nalign = 1

    # THIS STUDY ---------------------------------------------------------------
    ax.fill_between([1e-5, 1], [-0.5,-0.5], y2=[5.25, 5.25], color="#fff4e0", zorder=-10)
             
    # Plot actual values 
    for (d, e) in zip(dummyind_mine, plot_order_mine)
        x_thermal = [f_thermal[f_thermal.Exp.==e, :].Min[1], f_thermal[f_thermal.Exp.==e, :].Max[1]]
        x_both = [f_both[f_both.Exp.==e, :].Min[1], f_both[f_both.Exp.==e, :].Max[1]]
        y = fill!(similar(x_thermal), d)  # spaces them out in y axis

        # now actually plot them 
        ax.plot(x_thermal, y, linewidth=lw, color="white", zorder=4)#, solid_capstyle="round")
        ax.plot(x_thermal, y, linewidth=lw, color=colors[d], zorder=5, alpha=0.8)#, solid_capstyle="round")
        ax.scatter(x_thermal, y, color=scatcolors[d], zorder=6, marker=m, s=msize)
        ax.plot(x_both, y, linewidth=lw, color="white", zorder=4)
        ax.plot(x_both, y, linewidth=lw, color=colors[d], zorder=5, alpha=0.8)
        ax.scatter(x_both, y, color=scatcolors[d], zorder=6, marker=m, s=msize)
    end
    # mean cases
    ax.plot([f_mean_thermal, f_mean_thermal], [0, 4.5], color=myorange, zorder=-7)
    ax.scatter(f_mean_thermal, 0, color=myorange, zorder=5, marker=m, s=msize)
    ax.plot([f_mean_both, f_mean_both], [0, 4.5], color=myorange, zorder=-7) 
    ax.scatter(f_mean_both, 0, color=myorange, zorder=5, marker=m, s=msize)

    # text
    ax.text(1e-5, 5, "This study", color="black", fontsize=24, va="top")
    ax.text(2e-4, 4.75, "Thermal escape only")
    ax.text(6e-3, 4.75, "Thermal + non-thermal escape")
    ax.text(nalign, 0, "Standard atm.", color=myorange, ha=nameha)
    ax.text(f_mean_thermal-0.0003, 0, L"\overline{f}="*"$(round(f_mean_thermal, digits=3))", 
            color=myorange, ha="right", va="center")
    ax.text(f_mean_both-0.01, 0, L"\overline{f}="*"$(round(f_mean_both, digits=3))",
            color=myorange, ha="right", va="center")
    ax.text(nalign, 1, L"$\overline{T}_{surface} \pm 25$%", 
            color=colors[1], va="center", ha=nameha)
    ax.text(nalign, 2, L"$\overline{T}_{tropo} \pm 25$%", 
            color=colors[2], va="center", ha=nameha)
    ax.text(nalign, 3, L"$\overline{T}_{exobase} \pm 25$%", 
            color=colors[3], va="center", ha=nameha)
    ax.text(nalign, 4, L"Water 1-100 $\mathrm{\mu}$m", 
            color=colors[4], va="center", ha=nameha)

    # PAST STUDIES -------------------------------------------------------------
    for (d, e) in zip(dummyind_others, plot_order_others)
        x = toplot_others[e]
        y = fill!(similar(toplot_others[e]), d)
        c = colors[d]
        ax.scatter(x, y, linewidth=2, color=c, marker=m, s=msize, zorder=6)
    end

    ax.plot([1.15e-2, 6.8e-2], [6, 6], color=colors[6], linewidth=lw, linestyle=":", 
             alpha=0.4, zorder=4) 
    ax.plot(toplot_others["Krasnopolsky 2002"], fill!(similar(toplot_others["Krasnopolsky 2002"]), 7), 
            color=colors[7], linewidth=lw, alpha=0.4, zorder=4)
    ax.plot([1.02e-2, 2e-1], [8, 8], color=colors[8], linewidth=lw, linestyle=":", 
             alpha=0.4, zorder=4) 

    # text
    ax.text(1e-5, 9, "Past studies", color="black", fontsize=24, va="bottom")

    # Clarke+ 2019
    ax.text(nalign, 6, "Clarke+ 2019", color=colors[6], ha=nameha, va="center")
    ax.text(toplot_others["Clarke+ 2019"][1], 6.3, "Ls 280", color=colors[6], 
            ha="center", va="center", fontsize=14)
    ax.text(toplot_others["Clarke+ 2019"][2], 6.3, "Ls 318", color=colors[6], 
            ha="center", va="center", fontsize=14)
    ax.text(7e-2, 6, "?", color=colors[6], va="center", ha="center")
    ax.text(0.01, 6, "?", color=colors[6], va="center", ha="center")

    # Kras 2002
    ax.text(nalign, 7, "Kras. 2002", color=colors[7], ha=nameha, va="center")
    ax.text(toplot_others["Krasnopolsky 2002"][1], 7.3, "Solar min", 
            color=colors[7], ha="right", fontsize=14)
    ax.text(toplot_others["Krasnopolsky 2002"][2], 7.3, "mean", color=colors[7], 
            ha="center", fontsize=14)
    ax.text(toplot_others["Krasnopolsky 2002"][3], 7.3, "max", color=colors[7], 
            ha="center", fontsize=14)

    # Kras 2000
    ax.text(nalign, 8, "Kras. 2000", color=colors[8], ha=nameha, va="center")
    ax.text(toplot_others["Krasnopolsky 2000"][1], 8.3, "Diffusion profile 1", 
            color=colors[8], ha="center", fontsize=14)
    ax.text(toplot_others["Krasnopolsky 2000"][2], 8.3, "profile 2", 
            color=colors[8], ha="center", fontsize=14)
    ax.text(0.18, 8, "?", color=colors[8], va="center", ha="center")
    ax.text(1e-2, 8, "?", color=colors[8], va="center", ha="center")

    # Yung+ 1988
    ax.text(nalign, 9, "Yung\n1988", color=colors[9], ha=nameha, va="center")
    ax.text(0.22, 9, "Standard model (thermal escape only)", color=colors[9], ha="right", va="center", fontsize=14)
    
    # PLOT CONFIGURATION -------------------------------------------------------
    plt.margins(x=0.01, y=0.01)
    ax.set_xlabel(L"Fractionation factor $f$", fontsize=20)
    ax.set_xscale("log")
    ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax.set_yticks(collect(1:1:10))
    ax.tick_params(which="both", labelsize=18)
    ax.set_xlim([1e-5, 1])
    # save to the two useful folders
    savefig(base*"f_results.png", bbox_inches="tight")
    savefig(results_dir*"ALL STUDY PLOTS/f_results.png", bbox_inches="tight")
end

# Do my things =================================================================

# base = results_dir*"VarWaterTemp-HDO250/"
base = results_dir*"VarWaterTemp/"

plot_results_caltech_together(base)