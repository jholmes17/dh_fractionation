################################################################################
# plot_f_results_yung.jl
# TYPE: (3) Analysis - optional
# DESCRIPTION: Calculates f for a reproduction of the main case from Yung+1988.
#
# Eryn Cangi
# Created 31 May 2019
# Last edited: 20 July 2020
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


do_it_for_yung_case()
