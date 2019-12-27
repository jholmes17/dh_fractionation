################################################################################
# plot_water_profs_panel.jl
# TYPE: Supporting (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: plot three water profiles of high, mean, and low water in a
# 3-panel plot. For work on a poster in 2019 for various conferences including 
# the 1st MAVEN PSG of the year.
# 
# Eryn Cangi
# August 2019 (updated)
# Currently tested for Julia: 0.7
################################################################################

using PyPlot
using PyCall
using HDF5
using LaTeXStrings

# fundamental constants ========================================================
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg

# mars parameters
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm

const fullspecieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
                         :O1D, :H, :N2, :Ar, :CO2pl, :HOCO,
                         # species for deuterium chemistry:
                         :HDO, :OD, :HDO2, :D, :DO2, :HD, :DOCO];
specieslist=fullspecieslist;

# Functions ====================================================================

function get_ncurrent(readfile)
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

function n_tot(n_current, z)
    # get the total number density at a given altitude
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
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

function plot_3panel(Tp):
    #=
    Plots a 3-panel figure of water profiles, for low, nominal, and high water
    mixing ratios. Kind of deprecated in favor of plot_by_prum.

    Tp: a vector of 3 temperatures, floating point, at surface, tropopause, and 
        exobase.
    =#
    alt = (0:2e5:200e5)
    Temp(z::Float64) = Tpiecewise(z, Tp[1], Tp[2], Tp[3])

    # used in combination w/n_current. Gets index corresponding to a given altitude
    # DO NOT MOVE: has to be here because alt needs to be defined first
    n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])

    # CREATE THE WATER AND HDO PROFILES ============================================
    rf = "converged_standardwater.h5"#"/data/GoogleDrive/Phys/LASP/chaffincode-working/converged_standardwater.h5"
    n_current = get_ncurrent(rf)
    alt = h5read(rf,"n_current/alt")

    # calculate H2O saturation vapor pressure
    # first parenthetical converts whole expression to 1/cm^3, b/c 2nd parenthetical is in mmHg
    Psat_npv(T::Float64) = (133.3/(10^6 * boltzmannK * T))*(10^(-2445.5646/T + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 - 6.757169)) 
    H2Osat = map(x->Psat_npv(x), map(Temp, alt)) # get vertical saturation profile in 1/cm^3 for atmosphere
    H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP of water as a fraction of total atmo
    H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)] # slice from start to first minimum of saturation fraction
    H2Oinitfrac = [H2Oinitfrac; fill(minimum(H2Osatfrac), length(alt)-2-length(H2Oinitfrac))] # sets SVP to minimum for all altitudes above first altitude valued at minimum
    H2Oinitfrac1 = deepcopy(H2Oinitfrac)
    H2Oinitfrac2 = deepcopy(H2Oinitfrac)
    H2Oinitfrac3 = deepcopy(H2Oinitfrac)

    H2Oinitfrac1[findall(x->x<30e5, alt)] .= 1e-3 # set all values below 30 km to 0.0001
    H2Oinitfrac2[findall(x->x<50e5, alt)] .= 1e-4
    H2Oinitfrac3[findall(x->x<50e5, alt)] .= 1e-5 

    for i in [1:length(H2Oinitfrac1);] # ensures profile never exceeds saturation
        H2Oinitfrac1[i] = H2Oinitfrac1[i] < H2Osatfrac[i+1] ? H2Oinitfrac1[i] : H2Osatfrac[i+1]
        H2Oinitfrac2[i] = H2Oinitfrac2[i] < H2Osatfrac[i+1] ? H2Oinitfrac2[i] : H2Osatfrac[i+1]
        H2Oinitfrac3[i] = H2Oinitfrac3[i] < H2Osatfrac[i+1] ? H2Oinitfrac3[i] : H2Osatfrac[i+1]  
    end

    # plotting 
    # make plots pretty
    rcParams = PyDict(matplotlib["rcParams"])
    rcParams["font.sans-serif"] = ["Laksaman"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    fig, ax = subplots(1, 3, sharex=true, sharey=true, figsize=(10,4))
    subplots_adjust(wspace=0, bottom=0.25)

    # println(size(ax))
    ax[1][:set_title]("High water")
    ax[1][:semilogx](H2Oinitfrac1, alt[2:end-1]/1e5, label="1e-3")
    ax[1][:set_xlabel]("Water mixing ratio")
    ax[1][:set_ylabel]("Altitude (km)")
    ax[1][:set_xticks]([10.0^(-10), 10.0^(-8), 10.0^(-6), 10.0^(-4)])

    ax[2][:set_title]("Mean water")
    ax[2][:semilogx](H2Oinitfrac2, alt[2:end-1]/1e5, label="1e-4")
    ax[2][:set_xticks]([10.0^(-10), 10.0^(-8), 10.0^(-6), 10.0^(-4)])

    ax[3][:set_title]("Low water")
    ax[3][:semilogx](H2Oinitfrac3, alt[2:end-1]/1e5, label="1e-5")

    ax[3][:set_xticks]([10.0^(-10), 10.0^(-8), 10.0^(-6), 10.0^(-4)])


    savefig("../Results/VarWaterTemp/water_profs_3panel.png")
    show()
end

function plot_by_prum(Tp, hygropause, waters)
    alt = 0:2e5:250e5
    Temp(z::Float64) = Tpiecewise(z, Tp[1], Tp[2], Tp[3])
    n_current = get_ncurrent("converged.h5")

    # modify n_current with deuterated species profiles ============================
    n_current[:HDO] = n_current[:H2O] * DH
    n_current[:OD] = n_current[:OH] * DH
    n_current[:HDO2] = n_current[:H2O2] * DH
    n_current[:D] = n_current[:H] * DH
    n_current[:DO2] = n_current[:HO2] * DH
    n_current[:HD] = n_current[:H2] * DH
    n_current[:DOCO] = n_current[:HOCO] * DH
    # add the new Jrates --the values will get populated automatically =============
    n_current[:JHDOtoHpOD] = zeros(length(alt))
    n_current[:JHDOtoDpOH] = zeros(length(alt))
    n_current[:JHDO2toOHpOD] = zeros(length(alt))
    n_current[:JHDOtoHDpO1D] = zeros(length(alt)) 
    n_current[:JHDOtoHpDpO] = zeros(length(alt))
    n_current[:JODtoOpD] = zeros(length(alt))
    n_current[:JHDtoHpD] = zeros(length(alt))
    n_current[:JDO2toODpO] = zeros(length(alt))
    n_current[:JHDO2toDO2pH] = zeros(length(alt))
    n_current[:JHDO2toHO2pD] = zeros(length(alt))
    n_current[:JHDO2toHDOpO1D] =  zeros(length(alt))
    n_current[:JODtoO1DpD]  =  zeros(length(alt))


    H2Osat = map(x->Psat(x), map(Temp, alt)) # array in #/cm^3 by altitude
    H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP as fraction of total atmo

    fig = figure(figsize=(6,4))
    ax = gca()
    sevencols = ["#74a9cf", "#3690c0", "#0570b0","#034e7b"]#,"#a6bddb","#d0d1e6","#74a9cf"]
    prum = [25, 10, 1, 0.1]#[0.1, 1, 10, 25]# 50, 100, 150]
    j = 1
    for MR in waters

        # set H2O SVP fraction to minimum for all alts above first time min is reached
        H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
        H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
                       fill(minimum(H2Osatfrac), length(alt)-2-length(H2Oinitfrac))]

        #make profile constant in the lower atmosphere (well-mixed)
        H2Oinitfrac[findall(x->x<hygropause, alt)] .= MR
 
        # restrict the initial water fraction to be below the saturation vapor pressure curve
        for i in [1:length(H2Oinitfrac);]
            H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
        end
        
        linez = ["-", "-", "-", "-"]#, "-","-", ":"]
        if prum[j] in [0.1, 1, 10]
            L = "$(round(prum[j]; sigdigits=1)) pr "* L"$\mathrm{\mu}$m"
        else
            L = ""
        end
        semilogx(H2Oinitfrac, alt[2:end-1]./1e5, label=L, color=sevencols[j], linestyle=linez[j], linewidth=2)

        H2O_per_cc = sum([MR; H2Oinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
        H2Oprum = (H2O_per_cc * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
        j += 1
    end
    better_plots(ax)
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22
    
    text(10.0^(-5.6), -20, "0.1", color=sevencols[4], ha="right")
    text(10.0^(-4.6), -20, "1", color=sevencols[3], ha="right")
    text(10.0^(-2.8), -20, "10", color=sevencols[2], ha="right")
    text(10.0^(-2.6), -20, "25+", color=sevencols[1])
    text(10^(-10.5), -20, L"Column pr $\mathrm{\mu}$m:", fontsize=20, color="#555555")
    xlabel(L"H$_2$O Mixing ratio")
    ylabel("Altitude (km)")
    ylim(-25, 255)
    xlim(right=10.0^(-1.5))
    yticks([0,50,100,150,200,250])
    show()
    savefig("../Results/VarWaterTemp/250kmgrid/Plots/water_profiles.png", bbox_inches="tight")
end

# Make the stacked profile, showing how it changes with higher pr μm
plot_by_prum([216.0, 108.0, 211.0], 50e5, [2.18e-3, 8.1e-4, 2.1e-5, 1.48e-6])