################################################################################
# calc_f_converged.jl
# TYPE: MAIN (analysis)
# WHICH: Equilibrium experiments
# DESCRIPTION: easily calculate f for a given converged file. 
#
# Eryn Cangi
# 31 May 2019
# Currently tested for Julia: 0.7
################################################################################

using PyPlot
using PyCall
using HDF5
using Printf
using LaTeXStrings

# fundamental constants
boltzmannK = 1.38e-23;    # J/K
bigG = 6.67e-11;          # N m^2/kg^2
mH = 1.67e-27;            # kg

# mars parameters
marsM = 0.1075*5.972e24;  # kg
radiusM = 3396e5;         # cm

# altitude grid
alt = 0:2e5:250e5
zmax=alt[end]

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
    calculate_f_yung("/home/emc/GDrive-CU/Research/Results/Yung-With-Old-Water/Case 1/2019 redo for assurance/converged_yung_case1.h5", 1.2e8)
end

# Standard functions ===========================================================
function get_ncurrent(readfile)
    #=
    Retrieves the matrix of species concentrations by altitude from an HDF5
    file containing a converged atmosphere.
    =#
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    return n_current
end

global rcParams = PyCall.PyDict(matplotlib."rcParams")

function better_plots(axob)
    axob.set_facecolor("#ededed")
    axob.grid(zorder=0, color="white", which="both")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 16
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
end

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

function sbc(species, oflux, temps)
    #=
    Defines boundary conditions for a variable O flux and temperature set.
    species: symbol; the species to get the boundary condition.
    oflux: a number such as 1.2e8
    temps: an array of surface temp, tropopause temp, exobase temp. 
    =#

    Temp(z::Float64) = Tpiecewise(z, temps[1], temps[2], temps[3])
    H2Osat = map(x->Psat(x), map(Temp, alt))
    HDOsat = map(x->Psat_HDO(x), map(Temp, alt))
    H_effusion_velocity = effusion_velocity(Temp(zmax),1.0)
    H2_effusion_velocity = effusion_velocity(Temp(zmax),2.0)
    D_effusion_velocity = effusion_velocity(Temp(zmax),2.0)
    HD_effusion_velocity = effusion_velocity(Temp(zmax),3.0)

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

function get_all_H_flux(readfile, oflux, temps)
    #=
    Produces an array of H flux out of the top of the atmosphere, for a given
    O flux and temperature array.
    readfile: a converged simulation, type .h5
    oflux: a value such as 1.2e8 representing the O flux boundary condition
    temps: an array of surface temp, tropopause temp, exobase temp.  
    =#
    n_current = get_ncurrent(readfile)
    inds = Dict(158 => 1, 205 => 2, 264 => 3)  # convert the temp to an index
    i = inds[Int(temps[3])]

    # Nonthermal ecsape velocities for temperatures: T_exo = 158K, 205K, 264K. 
    # using ratios of thermal/nonthermal from Kras 2010
    # v_nt = Dict(:H => [30, 52, 92], :H2 => [5, 6, 7], :D => [12, 15, 18],
    #             :HD => [1, 6, 14])  #in cm/s
    v_nt = Dict(:H => [3.8, 49.3, 106.5], :H2 => [1, 5.04, 11.5], :D => [8.6, 15.5, 24.2],
                :HD => [0.19, 3.9, 8.5])  #in cm/s

    # Calculate the H flux: sum of ([H] [#/cm^3]) * (H flux rate [cm/s])
    # and (H2 population @ 200 km) * 2*(H2 flux rate) and 
    # (HD population @ 200 km) * (HD flux rate)
    Hfluxes = (n_current[:H][end]*sbc(:H, oflux, temps)[2,2]
                  + 2*n_current[:H2][end]*sbc(:H2, oflux, temps)[2,2]
                  + n_current[:HD][end]*sbc(:HD, oflux, temps)[2,2]
                  # NONTHERMAL
                  # )  # uncomment to turn off nonthermal
                  + n_current[:H][end]*v_nt[:H][i]
                  + 2*n_current[:H2][end]*v_nt[:H2][i]
                  + n_current[:HD][end]*v_nt[:HD][i])
    return Hfluxes
end

function get_all_D_flux(readfile, oflux, temps)
    #=
    Gets value of D flux out of the atmosphere in a converged file.
    =#
    n_current = get_ncurrent(readfile)
    inds = Dict(158 => 1, 205 => 2, 264 => 3)  # convert the temp to an index
    i = inds[Int(temps[3])]

    # Nonthermal ecsape velocities for temperatures: T_exo = 158K, 205K, 264K. 
    # using ratios of thermal/nonthermal from Kras 2010
    # v_nt = Dict(:H => [30, 52, 92], :H2 => [5, 6, 7], :D => [12, 15, 18],
    #             :HD => [1, 6, 14])  #nonthermal escape velocities in cm/s
    v_nt = Dict(:H => [3.8, 49.3, 106.5], :H2 => [1, 5.04, 11.5], :D => [8.6, 15.5, 24.2],
                :HD => [0.19, 3.9, 8.5])  #in cm/s

    # Calculate the D flux: sum of ([D] [#/cm^3]) * (D flux rate [cm/s])
    # and (HD population @ 200 km) * (HD flux rate)
    Dfluxes = (n_current[:D][end]*sbc(:D, oflux, temps)[2,2]
               + n_current[:HD][end]*sbc(:HD, oflux, temps)[2,2]
               # NONTHERMAL
               # )  # uncomment to turn off nonthermal
               + n_current[:D][end]*v_nt[:D][i]
               + n_current[:HD][end]*v_nt[:HD][i])
    return Dfluxes
end

function get_thermal_H_flux(readfile, oflux, temps)
    #=
    Produces an array of H flux out of the top of the atmosphere, for a given
    O flux and temperature array.
    readfile: a converged simulation, type .h5
    oflux: a value such as 1.2e8 representing the O flux boundary condition
    temps: an array of surface temp, tropopause temp, exobase temp.  
    =#
    n_current = get_ncurrent(readfile)
    
    # Calculate the H flux: sum of ([H] [#/cm^3]) * (H flux rate [cm/s])
    # and (H2 population @ 200 km) * 2*(H2 flux rate) and 
    # (HD population @ 200 km) * (HD flux rate)
    Hfluxes = (n_current[:H][end]*sbc(:H, oflux, temps)[2,2]
                  + 2*n_current[:H2][end]*sbc(:H2, oflux, temps)[2,2]
                  + n_current[:HD][end]*sbc(:HD, oflux, temps)[2,2])
    return Hfluxes
end

function get_thermal_D_flux(readfile, oflux, temps)
    #=
    Gets value of D flux out of the atmosphere in a converged file.
    =#
    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of ([D] [#/cm^3]) * (D flux rate [cm/s])
    # and (HD population @ 200 km) * (HD flux rate)
    Dfluxes = (n_current[:D][end]*sbc(:D, oflux, temps)[2,2]
               + n_current[:HD][end]*sbc(:HD, oflux, temps)[2,2])
    return Dfluxes
end

function effusion_velocity(Texo::Float64, m::Float64)
    #=
    Returns effusion velocity for a species.
    Texo: temperature of the exobase (upper boundary) in K
    m: mass of one molecule of species in amu
    =#
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+zmax))
    vth=sqrt(2*boltzmannK*Texo/(m*mH))
    v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)
    return v
end

function n_tot(n_current, z)
    specieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
                   :O1D, :H, :N2, :Ar, :CO2pl, :HOCO,
                   # species for deuterium chemistry:
                   :HDO, :OD, :HDO2, :D, :DO2, :HD, :DOCO];
    # get the total number density at a given altitude
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

Psat_HDO(T::Float64) = ((133.3*1e-6)/(boltzmannK * T))*(10^(-2445.5646/T 
                         + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 
                         - 6.757169))

Psat(T::Float64) = (133.3/(10^6 * boltzmannK * T))*(10^(-2445.5646/T 
                    + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 
                    - 6.757169))

function calculate_f(thefile, temps, Oflux, escape_type)
    #=
    A function to calculate f or a single simulation, if you just want to check 
    the value of f real quick.
    =#
    ncur = get_ncurrent(thefile)

    if escape_type=="thermal"
        Hf = get_thermal_H_flux(thefile, Oflux, temps)
        Df = get_thermal_D_flux(thefile, Oflux, temps)
    elseif escape_type=="both"
        Hf = get_all_H_flux(thefile, Oflux, temps)
        Df = get_all_D_flux(thefile, Oflux, temps)
    else
        println("Invalid escape type: $(escape_type)")
    end
    return 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1])
end

function search_subfolders(path, key)
    #=
    Searches the top level subfolders within path for all folders mathing a certain regex given by key.
    =#
    folders = []
    for (root, dirs, files) in walkdir(path)
        if root==path
            for dir in dirs
                push!(folders, joinpath(root, dir)) # path to directories
            end
        end
    end

    wfolders = filter(x->occursin(key, x), folders)
    return wfolders
end

# Main =========================================================================
function make_my_result_plot(escape_type)
    #=
    escape_type: "thermal" or "both". Makes a plot of fractionation factor by 
                 experiment for escape_tpe of atmospheric escape (thermal, or 
                 thermal + nonthermal).
    =#
    n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
    Oflux = 1.2e8# whatever Oflux is, in format "1.2e8" cm^-2 s^-1
    base = "/home/emc/GDrive-CU/Research/Results/VarWaterTemp/SolarMean/"

    meanTs = 216.0
    downTs = 162.0
    upTs = 270.0
    meanTt = 108.0
    downTt = 81.0
    upTt = 135.0
    meanTe = 205.0
    downTe = 158.0
    upTe = 264.0

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
                flist_mine[row, 3] = calculate_f(filetouse, t, Oflux, escape_type)
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

    savefig("../Results/ALL STUDY PLOTS/f-results-plot-$(escape_type).png", bbox_inches="tight")
end

function plot_results_caltech()
    #=
    escape_type: "thermal" or "both". Makes a plot of fractionation factor by 
                 experiment for escape_type of atmospheric escape (thermal, or 
                 thermal + nonthermal).
    =#
    n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
    Oflux = 1.2e8# whatever Oflux is, in format "1.2e8" cm^-2 s^-1
    base = "/home/emc/GDrive-CU/Research/Results/VarWaterTemp/SolarMean/"

    meanTs = 216.0
    lowTs = 162.0
    hiTs = 270.0
    meanTt = 108.0
    lowTt = 81.0
    hiTt = 135.0
    meanTe = 205.0
    lowTe = 158.0
    hiTe = 264.0

    # get high and low f in surfaces
    f_surf = Array{Float64}(undef, 3, 3) #array to store f
    f_surf[:, 1] = [lowTs, meanTs, hiTs]
    i = 1
    for Ts in f_surf[:, 1]
        # construct file name
        F = "temp_$(Int(Ts))_$(Int(meanTt))_$(Int(meanTe))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_surf[i, 2] = calculate_f(filetouse, [Ts, meanTt, meanTe], Oflux, "thermal")
        f_surf[i, 3] = calculate_f(filetouse, [Ts, meanTt, meanTe], Oflux, "both")
        i += 1
    end

    # get high and low f in tropopause
    f_tropo = Array{Float64}(undef, 3, 3) #array to store f
    f_tropo[:, 1] = [lowTt, meanTt, hiTt]
    i = 1
    for Tt in f_tropo[:, 1]
        # construct file name
        F = "temp_$(Int(meanTs))_$(Int(Tt))_$(Int(meanTe))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_tropo[i, 2] = calculate_f(filetouse, [meanTs, Tt, meanTe], Oflux, "thermal")
        f_tropo[i, 3] = calculate_f(filetouse, [meanTs, Tt, meanTe], Oflux, "both")
        i += 1
    end

    # get high and low f in exobase
    f_exo = Array{Float64}(undef, 3, 3) #array to store f
    f_exo[:, 1] = [lowTe, meanTe, hiTe]
    i = 1
    for Te in f_exo[:, 1]
        # construct file name
        F = "temp_$(Int(meanTs))_$(Int(meanTt))_$(Int(Te))"
        filetouse = base*F*"/converged_"*F*".h5"

        # calculate f
        f_exo[i, 2] = calculate_f(filetouse, [meanTs, meanTt, Te], Oflux, "thermal")
        f_exo[i, 3] = calculate_f(filetouse, [meanTs, meanTt, Te], Oflux, "both")
        i += 1
    end

    # get high and low f in water
    f_water = Array{Float64}(undef, 7, 3) #array to store f
    waterfolders = search_subfolders(base, r"water_[0-9]+\.[0-9]+e-[0-9]")
    i = 1

    for w in waterfolders

        # get the experiment name
        waterexp = match(r"water_[0-9]+\.[0-9]+e-[0-9]", w).match
        watermr = match(r"[0-9]+\.[0-9]+e-[0-9]", waterexp).match
        # construct file name
        filetouse = w * "/converged_"*waterexp*".h5"

        # calculate f
        f_water[i, 1] = parse(Float64, watermr)
        f_water[i, 2] = calculate_f(filetouse, [meanTs, meanTt, meanTe], Oflux, "thermal")
        f_water[i, 3] = calculate_f(filetouse, [meanTs, meanTt, meanTe], Oflux, "both")
        i += 1
    end

    skey = L"T$_{surface}$"
    tkey = L"T$_{tropopause}$"
    ekey = L"T$_{exobase}$"

    # iterate through all the detailed files and make an array of fractionation factors by experiment
    mn = base * "temp_216_108_205/converged_temp_216_108_205.h5"

    toplot_thermal = Dict("Mean conditions"=>[calculate_f(mn, [meanTs, meanTt, meanTe], Oflux, "thermal")],
                  ekey=>[minimum(f_exo[:, 2]), maximum(f_exo[:, 2])],
                  tkey=>[minimum(f_tropo[:, 2]), maximum(f_tropo[:, 2])],
                  skey=>[minimum(f_surf[:, 2]), maximum(f_surf[:, 2])],
                  "Water vapor"=>[minimum(f_water[:, 2]), maximum(f_water[:, 2])],
                  "Clarke+ 2019"=>[0.016, 0.047],
                  "Krasnopolsky 2002"=>[0.055, 0.082, 0.167],
                  "Krasnopolsky 2000"=>[0.016, 0.135],
                  "Yung+ 1988"=>[0.32]
             )

    toplot_both = Dict("Mean conditions"=>[calculate_f(mn, [meanTs, meanTt, meanTe], Oflux, "both")],
                  ekey=>[minimum(f_exo[:, 3]), maximum(f_exo[:, 3])],
                  tkey=>[minimum(f_tropo[:, 3]), maximum(f_tropo[:, 3])],
                  skey=>[minimum(f_surf[:, 3]), maximum(f_surf[:, 3])],
                  "Water vapor"=>[minimum(f_water[:, 3]), maximum(f_water[:, 3])],
                  "Clarke+ 2019"=>[0.016, 0.047],
                  "Krasnopolsky 2002"=>[0.055, 0.082, 0.167],
                  "Krasnopolsky 2000"=>[0.016, 0.135],
                  "Yung+ 1988"=>[0.32]
             )
    
    # THERMAL ESCAPE ONLY PLOT -------------------------------------------------
    fig, ax = subplots(figsize=(8,6))
    for side in ["top", "left", "right"]# "bottom", 
        ax.spines[side].set_visible(false)
    end
    ax.tick_params(labelleft=false, left=false)
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 16
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    # Which things to plot 
    plotorder = [skey, tkey, ekey, "Water vapor", "Yung+ 1988"]
    dummyind = [1, 2, 3, 4, 5]
    colors = ["xkcd:faded orange", "xkcd:faded orange", "xkcd:faded orange", 
              "xkcd:faded orange", "#485468"]

    # Plot actual values 
    for (d, e) in zip(dummyind, plotorder)
        x = toplot_thermal[e]
        y = fill!(similar(toplot_thermal[e]), d)
        c = colors[d]

        if d==5
            ax.scatter(x, y, linewidth=2, color=c, marker="o", s=80, zorder=6)
        else
            ax.plot(x, y, linewidth=15, color=c, marker="s", zorder=5)
        end
    end
    ax.plot([toplot_thermal["Mean conditions"], toplot_thermal["Mean conditions"]], 
               [0.5, 5], color="black", zorder=-500)  # mean condition line 

    # text
    ax.text(0.002, 1, L"$\overline{T}_{surface} \pm 25$%", color=colors[1], va="center")
    ax.text(0.003, 2, L"$\overline{T}_{tropopause} \pm 25$%", color=colors[2], va="center")
    ax.text(0.025, 3, L"$\overline{T}_{exobase} \pm 25$%", color=colors[3], va="center")
    ax.text(0.0015, 4, L"Water vapor 1-25 $\mu$m", color=colors[4], va="center")
    ax.text(0.26, 5, "Yung+ 1988", color=colors[5], ha="right", va="center")
    # mean value labels
    ax.text(0.0009 * (1-0.6), 4.8, "mean f = $(round(toplot_thermal["Mean conditions"][1], digits=3))\nthis study", va="top", ha="right")
    ax.plot([0.0004, 0.001], [4.4, 4.4], color="black")

    # basic plot things 
    plt.margins(y=0.05)
    ax.set_xlabel("Fractionation Factor", fontsize=20)
    ax.set_xscale("log")
    ax.set_yticks(collect(1:1:size(plotorder)[1]))
    ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax.set_yticklabels(plotorder)
    ax.tick_params(which="both", labelsize=18)  # for some reason this has to be here because editing rcParams in a function doesn't work if you do it more than once.
    plt.title("Thermal escape only")
    savefig("../Results/ALL STUDY PLOTS/f_results_plot_thermal_caltech.png", bbox_inches="tight")

    # NONTHERMAL + THERMAL ESCAPE PLOT -----------------------------------------

    # which things to plot
    fig, ax = subplots(figsize=(8,6))
    for side in ["top", "left", "right"]# "bottom", 
        ax.spines[side].set_visible(false)
    end
    ax.tick_params(labelleft=false, left=false)
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 16
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    plotorder = ["Clarke+ 2019", skey, tkey, ekey, "Water vapor", "Krasnopolsky 2002", 
                 "Krasnopolsky 2000"]
    dummyind = [1, 2, 3, 4, 5, 6, 7]
    colors = ["#3c6843", "xkcd:faded orange", "xkcd:faded orange", "xkcd:faded orange", 
              "xkcd:faded orange", "#682f51", "#99801d"]

    # plot values
    for (d, e) in zip(dummyind, plotorder)
        x = toplot_both[e]
        y = fill!(similar(toplot_both[e]), d)
        c = colors[d]

        if (d==1) || (d>5)
            ax.scatter(x, y, linewidth=2, color=c, marker="o", s=80, zorder=6)
        else
            ax.plot(x, y, linewidth=15, color=c, marker="s", zorder=5)
        end
    end
    ax.plot([toplot_both["Mean conditions"], toplot_both["Mean conditions"]], 
           [1.5, 5.5], color="black", zorder=-500)  # mean condition line 
    ax.plot(toplot_both["Krasnopolsky 2002"], fill!(similar(toplot_both["Krasnopolsky 2002"]), 6), 
            color=colors[6], linewidth=15, alpha=0.4, marker="o", zorder=4)
    
    # text
    ax.text(0.055, 1, "Clarke+ 2019, inferred from observations", color=colors[1], ha="left", va="center")
    ax.text(0.05, 2, L"$\overline{T}_{surface} \pm 25$%", color=colors[2], va="center")
    ax.text(0.08, 3, L"$\overline{T}_{tropopause} \pm 25$%", color=colors[3], va="center")
    ax.text(0.08, 4, L"$\overline{T}_{exobase} \pm 25$%", color=colors[4], va="center")
    ax.text(0.05, 5, L"Water vapor 1-25 $\mu$m", color=colors[5], va="center")
    ax.text(0.2, 6, "Krasnopolsky 2002", color=colors[6], ha="left", va="center")
    ax.text(0.18, 7, "Krasnopolsky 2000", color=colors[7], ha="left", va="center")
    # mean stuff
    ax.text(0.037 * (1-0.2), 6, "mean f = $(round(toplot_both["Mean conditions"][1], digits=3))\nthis study", va="top", ha="right")
    ax.plot([0.03, 0.038], [5.4, 5.4], color="black")
    ax.text(toplot_both["Clarke+ 2019"][1], 0.5, "Ls 280", color=colors[1], ha="left", va="center")
    ax.text(toplot_both["Clarke+ 2019"][2], 0.5, "Ls 318", color=colors[1], ha="left", va="center")
    ax.text(toplot_both["Krasnopolsky 2000"][1], 7.3, "Diffusion profile 1", color=colors[7], ha="center", fontsize=14)
    ax.text(toplot_both["Krasnopolsky 2000"][2], 7.3, "profile 2", color=colors[7], ha="center", fontsize=14)
    ax.text(toplot_both["Krasnopolsky 2002"][1], 6.3, "Solar min", color=colors[6], ha="right", fontsize=14)
    ax.text(toplot_both["Krasnopolsky 2002"][2], 6.3, "mean", color=colors[6], ha="center", fontsize=14)
    ax.text(toplot_both["Krasnopolsky 2002"][3], 6.3, "max", color=colors[6], ha="center", fontsize=14)

    
    plt.margins(y=0.15)
    ax.set_xlabel("Fractionation Factor", fontsize=20)
    ax.set_xscale("log")
    ax.set_xticks([1e-2, 1e-1, 1e0])
    ax.set_yticks(collect(1:1:size(plotorder)[1]))
    ax.set_yticklabels(plotorder)
    ax.tick_params(which="both", labelsize=18)
    plt.title("Thermal + nonthermal escape")
    savefig("../Results/ALL STUDY PLOTS/f_results_plot_both_caltech.png", bbox_inches="tight")
end

# Do my things =================================================================
plot_results_caltech()
# make_my_result_plot("thermal")
# make_my_result_plot("both")