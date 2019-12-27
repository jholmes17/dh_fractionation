################################################################################
# sep_atomic_molec_fluxes.jl
# TYPE: SUPPORTING (Troubleshooting?)
# WHICH: Equilibrium experiments
# DESCRIPTOIN: An annoyingly long file to calculate the separate H, H2, D, HD 
# fluxes in a converged atmosphere file.
#
# Eryn Cangi 
# April 2019
# Currently tested for Julia: 0.7
################################################################################

using PyPlot
using HDF5

# CONSTANTS ====================================================================

# fundamental constants
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg

# mars parameters
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm

DH = 5.5 * 1.6e-4        # SMOW value from Yung 1988

alt = (0:2e5:200e5)
zmax = alt[end]

# FUNCTIONS ====================================================================
function get_ncurrent(readfile)
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
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

Temp(z::Float64) = Tpiecewise(z, 192.0, 110.0, 199.0)

function n_tot(n_current, z)
    # get the total number density at a given altitude
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

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

function speciesbcs(species)
    get(speciesbclist,
        species,
        ["f" 0.; "f" 0.])
end

function sep_atomic_molec_fluxes(readfile)
    #=
    Produces an array of H and D fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
    =#
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    thematrix = mydata["n_current"]["n_current_mat"]
    close(mydset)

    

    n_current = get_ncurrent(readfile)

    Hfluxes = n_current[:H][end]*speciesbcs(:H)[2,2]
    H2fluxes = n_current[:H2][end]*speciesbcs(:H2)[2,2]
    Dfluxes = n_current[:HD][end]*speciesbcs(:HD)[2,2]
    HDfluxes = n_current[:D][end]*speciesbcs(:D)[2,2]

    return Hfluxes, H2fluxes, Dfluxes, HDfluxes
end

# FILENAME =====================================================================
experimentdir = "/home/emc/GDrive-CU/Research/Results/VarWaterTemp/"  
FNext = "temp_192_110_199"
filename = experimentdir*FNext*"/converged_standardwater_D_"*FNext*".h5"

n_current = get_ncurrent(filename)
println("H pop: ", n_current[:H][end])
println("H2 pop: ", n_current[:H2][end])
println("D pop: ", n_current[:D][end])
println("HD pop: ", n_current[:HD][end])
println()
n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])

specieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO, :O1D, :H, 
               :N2, :Ar, :CO2pl, :HOCO, :HDO, :OD, :HDO2, :D, :DO2, :HD, :DOCO];

# CODE NECESSARY FOR BOUNDARY CONDITIONS =======================================
Psat(T::Float64) = ((133.3*1e-6)/(boltzmannK * T))*(10^(-2445.5646/T + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 - 6.757169))
H2Osat = map(x->Psat(x), map(Temp, alt)) # array in #/cm^3 by altitude
H2Osatfrac = H2Osat./map(z->n_tot(n_current, z),alt)  # get SVP as fraction of total atmo
# set H2O SVP fraction to minimum for all alts above first time min is reached
H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)),H2Osatfrac), 0)]
H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
               fill(minimum(H2Osatfrac), length(alt)-2-length(H2Oinitfrac))]
               
HDOsat = map(x->Psat_HDO(x), map(Temp, alt))

H_effusion_velocity = effusion_velocity(Temp(zmax),1.0)
H2_effusion_velocity = effusion_velocity(Temp(zmax),2.0)
D_effusion_velocity = effusion_velocity(Temp(zmax),2.0)
HD_effusion_velocity = effusion_velocity(Temp(zmax),3.0)

println("H v_eff: ", H_effusion_velocity)
println("H2 v_eff: ", H2_effusion_velocity)
println("D v_eff: ", D_effusion_velocity)
println("HD v_eff: ", HD_effusion_velocity)
println()

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

Hf, H2f, Df, HDf = sep_atomic_molec_fluxes(filename)
println("H flux: ", Hf)
println("H2 flux: ", H2f)
println("D flux: ", Df)
println("HD flux: ", HDf)
