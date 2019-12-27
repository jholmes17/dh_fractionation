################################################################################
# reproductions_plot.jl
# TYPE: SUPPORTING (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: Makes a plot of reproductions of past studies. Is a separate 
# file for simplicity and because it allows more easy access to modify the way
# I calculate the reproduction of Krasnopolsky 2002 (which can be done either
# as the models were run, or by including his nonthermal escape velocities)
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
    # rcParams["font.sans-serif"] = ["Louis George Caf?"]
    # rcParams["font.monospace"] = ["FreeMono"]
    # rcParams["font.size"] = 16
    # rcParams["axes.labelsize"]= 20
    # rcParams["xtick.labelsize"] = 18
    # rcParams["ytick.labelsize"] = 18
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
    ztropo = 90e5  # height of the tropopause top
    
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
        return Texo - (Texo - Ttropo)*exp(-((z-ztropo)^2)/(11.4e10*Texo))
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

function get_H_flux(readfile, oflux, temps)
    #=
    Produces an array of H flux out of the top of the atmosphere, for a given
    O flux and temperature array.
    readfile: a converged simulation, type .h5
    oflux: a value such as 1.2e8 representing the O flux boundary condition
    temps: an array of surface temp, tropopause temp, exobase temp.  
    =#

    inds = Dict(200 => 1, 270 => 2, 350 => 3)  # convert the temp to an index
    i = inds[temps[3]]

    # nonthermal escape velocities for Krasnopolsky 2002 reproduction.
    # Each species has a value recorded at T = 200K, 270K, and 350K.
    v_nt = Dict(:H => [38, 56, 89], :H2 => [12.9, 18.2, 28], :D => [17, 24, 37],
                :HD => [8.2, 11.5, 17.7])  #nonthermal escape velocities in cm/s

    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (H population @ 200 km) * (D flux rate)
    # and (H2 population @ 200 km) * 2*(H2 flux rate) and 
    # (HD population @ 200 km) * (HD flux rate)
    # INCLUDES NONTHERMAL!!!
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

function get_D_flux(readfile, oflux, temps)
    #=
    Gets value of D flux out of the atmosphere in a converged file.
    =#

    inds = Dict(200 => 1, 270 => 2, 350 => 3)  # convert the temp to an index
    i = inds[temps[3]]

    v_nt = Dict(:H => [38, 56, 89], :H2 => [12.9, 18.2, 28], :D => [17, 24, 37],
                :HD => [8.2, 11.5, 17.7])  #nonthermal escape velocities in cm/s

    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (D population @ 200 km) * (D flux rate)
    # and (HD population @ 200 km) * (HD flux rate)
    # INCLUDES NONTHERMAL!!!
    Dfluxes = (n_current[:D][end]*sbc(:D, oflux, temps)[2,2]
               + n_current[:HD][end]*sbc(:HD, oflux, temps)[2,2]
               # NONTHERMAL
               # )  # uncomment to turn off nonthermal
               + n_current[:D][end]*v_nt[:D][i]
               + n_current[:HD][end]*v_nt[:HD][i])
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

function calculate_f(thefile, temps, Oflux)
    #=
    A function to calculate f or a single simulation, if you just want to check 
    the value of f real quick.
    =#
    ncur = get_ncurrent(thefile)
    Hf = get_H_flux(thefile, Oflux, temps)
    Df = get_D_flux(thefile, Oflux, temps)
    return 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1])
end

# Main =========================================================================

function make_reproduction_plot(fmin, fmean, fmax, therm)
    #=
    fmin: f for the solar minimum reproduction of Krasnopolsky 2002.
    fmean: f mean for the solar mean of same
    fmax: you know the drill
    therm: whether to reproduce his results using just thermal escape, or by 
           also accounting for nonthermal escape (in a ham-fisted way, by just
           multiplying his nonthermal escape velocities by the species 
           concentration)
    =#
    n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
    Oflux = 1.2e8 

    # REPRODUCTION PLOT ========================================================
    flist_repro = [["past", "Replication", 0.26],
                   ["past", "Yung+1988", 0.32],
                   ["past", "Replication", round(fmin, sigdigits=2)],
                   ["past", L"Kras. 2002 $\odot$ min", 0.055],
                   ["past", "Replication", round(fmean, sigdigits=2)],
                   ["past", L"Kras. 2002 $\odot$ mean", 0.082],
                   ["past", "Replication", round(fmax, sigdigits=2)],
                   ["past", L"Kras. 2002 $\odot$ max", 0.167]]

    flist_repro = permutedims(reshape(hcat(flist_repro...), (length(flist_repro[1]), length(flist_repro))))

    carray_repros = ["#74c476", "#74c476", "#a3bfdc", "#a3bfdc", "#8c96c6", 
                     "#8c96c6", "#8c6bb1", "#8c6bb1"]

    fig, ax = subplots(figsize=(7,7))
    better_plots(ax)
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 16
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    ax.barh(collect(1:1:size(flist_repro)[1]), flist_repro[:, 3], 0.9, color=carray_repros, zorder=10)
    ax.set_xscale("log")
    ax.set_xlabel("Fractionation Factor", fontsize=20)
    ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.set_yticks(collect(1:1:size(flist_repro)[1]))
    ax.set_yticklabels(flist_repro[:, 2], fontsize=18)
    ax.tick_params(labelsize=18)
    for (i, v) in enumerate(flist_repro[:, 3])
        if v <= 2e-5
            m = 4
        else
            m = -0.05
        end

        printvalstyle = "$(v)"
        ax.text(v+m*v, i-0.1, printvalstyle, color="black", zorder=16,
                ha="right")
    end

    if therm=="thermal"
        ax.set_title("Study replications, thermal escape only", fontsize=22)
        savefig("../Results/ALL STUDY PLOTS/f-reproductions-plot-thermal.png", bbox_inches="tight")
    elseif therm=="nonthermal"
        ax.set_title("Study replications, incl. nonthermal escape", fontsize=22)
        savefig("../Results/ALL STUDY PLOTS/f-reproductions-plot-nonthermtoo.png", bbox_inches="tight")
    end
end

# Do my things =================================================================

filemin = "/home/emc/GDrive-CU/Research/Results/Kras2002-rep/temp_213_125_200/converged_temp_213_125_200.h5"
filemean = "/home/emc/GDrive-CU/Research/Results/Kras2002-rep/temp_213_125_270/converged_temp_213_125_270.h5"
filemax = "/home/emc/GDrive-CU/Research/Results/Kras2002-rep/temp_213_125_350/converged_temp_213_125_350.h5"
f_k2002min = calculate_f(filemin, [213.0, 125.0, 200.0], 1.2e8)
f_k2002mean = calculate_f(filemean, [213.0, 125.0, 270.0], 1.2e8)
f_k2002max = calculate_f(filemax, [213.0, 125.0, 350.0], 1.2e8)

println(f_k2002min)
println(f_k2002mean)
println(f_k2002max)

make_reproduction_plot(f_k2002min, f_k2002mean, f_k2002max)