################################################################################
# check_eq.jl
# TYPE: Supporting (Verification)
# WHICH: Equilibrium and perturbation experiments
# DESCRIPTION: Check if the atmosphere is in equilibrium by looking for
# Φ_O = 0.5Φ_{H+D}.
#
# Eryn Cangi
# 11 September 2018
# Currently tested for Julia: 0.7
################################################################################
using HDF5

# fundamental constants ========================================================
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm
DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Experiment type and loading the file =========================================
# get arguments to use for accessing files 
argarray = Any[ARGS[i] for i in 1:1:length(ARGS)] 
if argarray[1] == "temp"
    FNext = "temp_$(argarray[2])_$(argarray[3])_$(argarray[4])"
else
    FNext = "$(argarray[1])_$(argarray[2])"
end

# load readfile and altitude array
# lead = "/data/GDrive-CU/Research/Results/"
lead = "/home/emc/GDrive-CU/Research/Results/"
rf = lead * FNext * "/converged_standardwater_D_" * FNext * ".h5"
println(rf)
alt = h5read(rf,"n_current/alt")

const zmax = alt[end];

# Functions ====================================================================
function speciesbcs(species)
    get(speciesbclist,
        species,
        ["f" 0.; "f" 0.])
end

function Tpiecewise(z::Float64, Tsurf, Ttropo, Tinf, lapserate=-1.4e-5, ztropo=90e5)
    #=
    DO NOT MODIFY! If you want to change the temperature, define a
    new function or select different arguments and pass to Temp(z)

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 temperatures for altitudes
    >htropo=90km, fixed at Ttropo=125K between htropo and
    htropo-htropowidth=60km, and rising at a constant lapse rate
    (1.4K/km) below.

    z: an altitude in cm.
    returns: temperature in K at the requested altitude.
    =#
    # allow varying troposphere width
    ztropowidth = ztropo - (1/lapserate)*(Ttropo-Tsurf)
    if ztropowidth < 0   # in case a profile is nonphysical, let lapse rate vary
        ztropowidth = 30e5
        m = (ztropo/1e5 - ztropowidth/1e5) / (Ttropo - Tsurf)
        lapserate = (1/m)*1e-5
    end

    if z >= ztropo
        return Tinf - (Tinf - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Tinf))
    end
    if ztropo > z >= ztropo - ztropowidth
        return Ttropo
    end
    if ztropo-ztropowidth > z
        return Ttropo-lapserate*(ztropo-ztropowidth-z)
    end
end

function get_ncurrent(readfile)
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

function get_H_and_D_fluxes(readfile)
    #=
    Produces an array of H and D fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
    =#

    n_current = get_ncurrent(readfile)

    # Calculate the total H+D flux: sum of (H population @ 200 km) * (H flux rate)
    # and (H_2 population @ 200 km) * (H_2 flux rate), and deuterated species
    HDfluxes = (n_current[:H][end]*speciesbcs(:H)[2,2]       # H
                  + 2*n_current[:H2][end]*speciesbcs(:H2)[2,2]  # H
                  + n_current[:HD][end]*speciesbcs(:HD)[2,2]    # H
                  + n_current[:D][end]*speciesbcs(:D)[2,2]      # D
                  + n_current[:HD][end]*speciesbcs(:HD)[2,2]) # D
    return HDfluxes
end

function get_H_fluxes(readfile)
    #=
    Produces an array of H and D fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
    =#

    n_current = get_ncurrent(readfile)

    # Calculate the total H+D flux: sum of (H population @ 200 km) * (H flux rate)
    # and (H_2 population @ 200 km) * (H_2 flux rate), and deuterated species
    Hfluxes = (n_current[:H][end]*speciesbcs(:H)[2,2]
                  + 2*n_current[:H2][end]*speciesbcs(:H2)[2,2])
    return Hfluxes
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

Psat(T::Float64) = (133.3/(10^6 * boltzmannK * T))*(10^(-2445.5646/T + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 - 6.757169))

function LHDO(T)
    #=
    Latent heat of vaporization of HDO as a function of temperature in K.
    This analytical function was determined by fitting data from
    https://pubs.acs.org/doi/pdf/10.1021/cr60292a004 (Jancso 1974, page 734)
    and extrapolating. It is probably not accurate outside the range of the data,
    which was 0-100, but it shouldn't be too far off.

    The data was in cal/mol versus Celsius. We convert to kJ/mol below.
    Fit was done in Python in a separate script. parameters below are the output
    from said fit.
    =#
    a = -0.02806171415983739
    b = 89.51209910268079
    c = 11918.608639939 # this is cal/mol
    return (a*(T-b)^2 + c) / 239 # returns in kJ/mol
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

# Set up non-function code that is necessary ===================================
if argarray[1]=="temp"
    Temp(z::Float64) = Tpiecewise(z, parse(Float64,argarray[2]), 
                                     parse(Float64,argarray[3]), 
                                     parse(Float64,argarray[4]))
elseif argarray[1] == "orig"
    Temp(z::Float64) = Tpiecewise(z, 209.0, 125.0, 240.0)
else
    Temp(z::Float64) = Tpiecewise(z, 192.0, 110.0, 199.0)
end

# H2Osat is an array of H2O SVP in 1/cm^3, a number per volume, by altitude
H2Osat = map(x->Psat(x), map(Temp, alt))
HDOsat = map(x->Psat_HDO(x), map(Temp, alt))

H_effusion_velocity = effusion_velocity(Temp(zmax),1.0)
H2_effusion_velocity = effusion_velocity(Temp(zmax),2.0)
D_effusion_velocity = effusion_velocity(Temp(zmax),2.0)
HD_effusion_velocity = effusion_velocity(Temp(zmax),3.0)

const speciesbclist=Dict(
                :CO2=>["n" 2.1e17; "f" 0.],
                :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                :HDO=>["n" HDOsat[1]; "f" 0.],
                :O=>["f" 0.; "f" 1.2e8],
                :H2=>["f" 0.; "v" H2_effusion_velocity],
                :HD=>["f" 0.; "v" HD_effusion_velocity],
                :H=>["f" 0.; "v" H_effusion_velocity],
                :D=>["f" 0.; "v" D_effusion_velocity],
               );

# Calculate the flux ratio =====================================================
Of = 1.2e8
Hf = get_H_fluxes(rf)
HDf = get_H_and_D_fluxes(rf)

println("O flux: ", Of)
println("H+D flux: ", HDf)
println("H flux: ", Hf)
println("Ratio: ", HDf/Of)
println()
