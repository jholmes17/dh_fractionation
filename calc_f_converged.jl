################################################################################
# calc_f_converged.jl
# TYPE: MAIN (analysis)
# WHICH: Equilibrium experiments
# DESCRIPTION: easily calculate f for a given converged file. 
#
# Eryn Cangi
# 31 May 2019
# Currently tested for Julia: 0.7
# HAS NOT been tested to be working, though. Originally was run in Jupyter.
################################################################################

using HDF5

# fundamental constants
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg

# mars parameters
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm

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

function speciesbcs_VARIABLE(species, oflux, temps)
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

    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (H population @ 200 km) * (D flux rate)
    # and (H2 population @ 200 km) * 2*(H2 flux rate) and 
    # (HD population @ 200 km) * (HD flux rate)
    Hfluxes = (n_current[:H][end]*speciesbcs_VARIABLE(:H, oflux, temps)[2,2]
                  + 2*n_current[:H2][end]*speciesbcs_VARIABLE(:H2, oflux, temps)[2,2]
                  + n_current[:HD][end]*speciesbcs_VARIABLE(:HD, oflux, temps)[2,2])
    return Hfluxes
end

function get_D_flux(readfile, oflux, temps)
    #=
    Gets value of D flux out of the atmosphere in a converged file.
    =#
    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (D population @ 200 km) * (D flux rate)
    # and (HD population @ 200 km) * (HD flux rate)
    Dfluxes = (n_current[:D][end]*speciesbcs_VARIABLE(:D, oflux, temps)[2,2]
                  + n_current[:HD][end]*speciesbcs_VARIABLE(:HD, oflux, temps)[2,2])
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

function LHDO(T)
    #=
    Latent heat of vaporization of HDO as a function of temperature in K.
    Using Jancso 1974 (https://pubs.acs.org/doi/pdf/10.1021/cr60292a004)
    and Jacobson's Fundamentals of Atmospheric Modeling.

    How I determined this analytical function:
    1. Used Jancso 1974, pg. 734, to get data on the difference of latent heat 
       of evaporation (L_e) for HDO and H2O from 0 to 100°C. 
    2. Used equation 2.54 in Jacobson to get analytical values for L_e of H2O in 
       the same range (0-100°C)
    3. Combined the L_e of H2O with data from Jancso 1974 to get L_HDO.
    4. Noting that L_e expression is roughly linear, fit a line to the L_HDO
       values and extrapolated to Mars temperatures of interest (100 to 350K)

    The data from Jancso 1974 are in cal/mol vs Celsius. We convert to kJ/mol 
    below. 1 kJ = 239 cal => 0.004184 kJ/1 cal. Fit was done in Python in a 
    separate script. parameters below are the output from said fit.
    =#
      
    M = -11.368481383727467  
    a = -2.1942156994200985
    B = 11007.47142671782 # this is cal/mol
    return (M*(T-a) +B) * (0.004184 / 1) # returns in kJ/mol
end

function Psat_HDO(T)
    #=
    Analytical expression for saturation vapor pressure of HDO, using analytical
    latent heat for HDO determined from fitting to data from Jancso 1974.
    The boiling temperature for HDO (100.7°C, 373.85K) and standard pressure are
    used as reference temp/pressure. The gas constant for HDO was calculated as
    R_HDO = kB/m, which probably came from Pierrehumbert's text.

    returns the pressure in #/cm^3. The conversion from N/m^2 to #/cm^3 is 
    divide by N*m (that is, kT) and multiply by the conversion to cm from m 
    (1e-6). 

    Input
        T: a single temperature in K
    Output:
        Pressure in #/cm^3
   =#
   R_HDO = 434.8 # J/kgK
   L_HDO = LHDO(T) * 1000 / 0.019 # kJ/mol * J/kJ * mol / kg = J / kg
   T0 = 373.85 # boiling temp for liquid HDO
   P0 = 101325 # standard pressure
   Psat_pascals = P0*exp(-(L_HDO/R_HDO)*(1/T - 1/T0))
   # To convert to #/cm^3 from Pa:
   # Pa = Nm^-2 (1/(Nm)) * (1e-6 m^3 / 1 cm^3)
   return Psat_pascals*(1e-6)/(boltzmannK*T) # return in #/cm^3
end

Psat(T::Float64) = (133.3/(10^6 * boltzmannK * T))*(10^(-2445.5646/T 
                    + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 
                    - 6.757169))

alt = 0:2e5:200e5
zmax=alt[end]

n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])

function calculate_f(thefile, temps, Oflux)
    
    ncur = get_ncurrent(thefile)
    Hf = get_H_flux(thefile, Oflux, temps)
    Df = get_D_flux(thefile, Oflux, temps)
    println("file: ", thefile)
    println("f = ", 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
    println()
end

base = "/home/emc/GDrive-CU/Research/Results/"

# get arguments
args = Any[ARGS[i] for i in 1:1:length(ARGS)]
if args[1]=="temp"
    FNext = "temp_$(args[2])_$(args[3])_$(args[4])"
    for i in 2:1:length(args)
        args[i] = parse(Float64, args[i])
    end
    t = [args[2], args[3], args[4]]
elseif args[1]=="water"
    FNext = "water_$(args[2])"
    t = [190.0, 110.0, 200.0]
elseif args[1]=="dh"
    FNext = "dh_$(args[2])"
elseif args[1]=="Oflux"
    FNext = "Oflux_$(args[2])"
end


filetouse = base*"VarWaterTemp/"*FNext*"/converged_standardwater_D_"*FNext*".h5"
Oflux = 1.2e8# whatever Oflux is, in format "1.2e8" cm^-2 s^-1
calculate_f(filetouse, t, Oflux)