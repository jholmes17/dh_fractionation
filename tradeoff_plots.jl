################################################################################
# make_tradeoff_plots.jl
# TYPE: MAIN (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: make trade off plots showing H and D flux and other variables as 
# a function of random things such as water column. 

# Eryn Cangi
# 5 April 2019
# Currently tested for Julia: 0.7
################################################################################
using PyPlot
using HDF5
using LaTeXStrings
using PyCall
using PlotUtils
using JLD

# fundamental constants ========================================================
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm
DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Basic altitude stuff
alt = (0:2e5:250e5)
n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
const zmax = alt[end];

# Functions ====================================================================
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

function get_ncurrent(readfile)
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

function get_H_fluxes(readfile, oflux, temps)
    #=
    Produces an array of H and D fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
    =#

    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (H population @ 200 km) * (D flux rate)
    # and (H2 population @ 200 km) * 2*(H2 flux rate) and 
    # (HD population @ 200 km) * (HD flux rate)
    Hfluxes = (n_current[:H][end]*speciesbcs_var(:H, oflux, temps)[2,2]
                  + 2*n_current[:H2][end]*speciesbcs_var(:H2, oflux, temps)[2,2]
                  + n_current[:HD][end]*speciesbcs_var(:HD, oflux, temps)[2,2])
    return Hfluxes
end

function get_D_fluxes(readfile, oflux, temps)
    #=
    Gets value of D flux out of the atmosphere in a converged file.
    =#
    n_current = get_ncurrent(readfile)

    # Calculate the D flux: sum of (D population @ 200 km) * (D flux rate)
    # and (HD population @ 200 km) * (HD flux rate)
    Dfluxes = (n_current[:D][end]*speciesbcs_var(:D, oflux, temps)[2,2]
                  + n_current[:HD][end]*speciesbcs_var(:HD, oflux, temps)[2,2])
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

Psat(T::Float64) = (133.3/(10^6 * boltzmannK * T))*(10^(-2445.5646/T 
                    + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 
                    - 6.757169))

Psat_HDO(T::Float64) = (133.3/(10^6 * boltzmannK * T))*(10^(-2445.5646/T 
                    + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 
                    - 6.757169))

function speciesbcs_var(species, oflux, temps)

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

function better_plot_bg(axob)
    axob.set_facecolor("#ededed")
    axob.grid(zorder=0, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
end

# Special functions for this file ==============================================

searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function normalize(arr, base_i)
    normed_arr = arr ./ arr[base_i]
    return normed_arr
end

function normalize_val(arr, normval)
    normed_arr = arr ./ normval
    return normed_arr
end

function get_colors(L, cmap)
    #=
    Generates some colors based on a color map for use in plotting a bunch of
    lines all at once.
    L: number of colors to generate.
    cmap: color map name
    =#

    colors_dumb = [cgrad(Symbol(cmap))[x] for x in range(0, stop=1, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = red(colors_dumb[i])
        c[i, 2] = green(colors_dumb[i])
        c[i, 3] = blue(colors_dumb[i])
    end
    return c
end

function areadensity_to_micron_atm(numpercm2)
    # #/cm^2 * (cm^2/m^2) * (10μm-am/2.687e20 #/m^2)
    return numpercm2 * (1e4/1) * (10/2.687e20)
end

# Special plotting functions (apply to all main plotting funcs) ----------------

function DH_alt_prof_plot(DHproflist, exps, v, s, optext="", optlegend="")
    #=
    DHproflist: An array with D/H profiles in each row
    exps: a list of experiment identifiers as strings
    v: water, Oflux, or temp (for putting in correct folder)
    s: a subfolder "abs/" or "mr/" for either absolute abundances or mixing ratio
    optext: an optional extension that will be appended to v, for 
            the temperature case where it goes in the temp_plots 
            folder but the files also specify exobase, tropopause, etc.
    optlegend: an optional string to put into the legend
    =#

    # do the DH altitudinal profile plot
    # set up plot
    fig, ax = subplots(figsize=(6,4))
    better_plot_bg(ax)
    subplots_adjust(wspace=0, bottom=0.15)
    ax.set_xlabel("D/H ratio (in atomic D, H)")
    ax.set_ylabel("Altitude (km)")
    ax.set_yticks(ticks=collect(0:50:200))
    
    # generate colors
    c = get_colors(length(exps), "viridis")

    # do actual plotting
    if typeof(exps[1]) != String
        exps = [string(x) for x in exps]
    end
    for i in range(1, length=length(exps))
        ax.plot(DHproflist[i, :], alt[2:end-1]./1e5, zorder=10, color=c[i, :], 
                linewidth=3, label=optlegend*"="*exps[i])
    end

    # set savepath
    plotpath = lead*v*"_plots/"*s
    savepath = plotpath*v*optext*"_DH_prof.png"
    #legend(fontsize=12, bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
    savefig(savepath, bbox_inches="tight")

    # save it again but with log x-axis 
    xscale("log")
    xticks(rotation=45)
    savepath = plotpath*"/"*v*optext*"_DH_prof_LOG.png"
    savefig(savepath, bbox_inches="tight")
    close(fig)
end

function CO_O2_plot(xvals, ydict, xlab, pathkey, meanX, s, tempkey="")
    #=
    xvals: the variations within the experiment. list of strings
    ydict: the dictionary containing the data
    xlab: label for the x axis
    pathkey: "water", "Oflux" or "temp"
    meanX: the nominal value to plot
    s: a subfolder "abs/" for absolute abundances or "mr/" for mixing ratio
    tempkey: "_surface", "_tropopause", "_exobase"
    =#

    # Pretty ugly, but data structs are always better than if/else, right?
    paststudies = Dict("water"=>Dict("yung"=>15, "nair"=>[3,8.8]), # these are in ppm
                            "Oflux"=>Dict("yung"=>5, "nair"=>9), # indices
                            "temp_surface"=>Dict("yung"=>220, "nair"=>214), 
                            "temp_tropopause"=>Dict("yung"=>140, "nair"=>140), 
                            "temp_exobase"=>Dict("yung"=>364, "nair"=>288))

    lookupkey = pathkey*tempkey
    ystr = s == "abs/" ? L"Abundance, CO and O$_2$" : L"Mixing ratio, CO and O$_2$"
    
    # CO, O2, CO/O2 plot
    # set up plot
    fig, ax = subplots(figsize=(6,4))
    xlabel(xlab)
    if pathkey=="water"
        xscale("log")
    end
    xticks(rotation=45)
    medgray = "#444444"

    # CO and O2 mixin ratio axis
    ax.plot(xvals, ydict["CO"], marker="o", zorder=10, color="navy", label="CO", alpha=0.5) 
    ax.plot(xvals, ydict["O2"], marker="D", zorder=10, color="navy", label=L"O$_2$", alpha=0.5)
    
    ax.axvline(meanX, color=medgray, zorder=10)
    nomtexty = Dict("temp_exobase"=>8.5e-4, "temp_surface"=>2e-3, 
                    "temp_tropopause"=>1.2e-3, "water"=>4e-3, "Oflux"=>6e-4)
    ax.text(meanX, nomtexty[lookupkey], "Nomimal\ncase", color=medgray)
    ax.set_facecolor("#ededed")
    ax.tick_params(axis="y", labelcolor="navy", which="both")
    ax.set_ylabel(ystr, color="navy")
    ax.set_yscale("log")
    # ax.legend()
    for side in ["top", "bottom", "left", "right"]
        ax.spines[side].set_visible(false)
    end
    
    # CO/O2 axis
    ax2 = ax.twinx() # set up the other axis
    ax2.axhline(0.6, color="red", linestyle="--")
    obstxtx = Dict("temp_exobase"=>170, "temp_tropopause"=>80, "temp_surface"=>250, 
                   "water"=>1e-1, "Oflux"=>0)
    obstxty = Dict("temp_exobase"=>0.55, "temp_tropopause"=>0.7, "temp_surface"=>0.55, 
                   "water"=>0.55, "Oflux"=>0.7)
    ax2.text(obstxtx[lookupkey], obstxty[lookupkey], "Obs.", color="red", va="top")
    ax2.grid(zorder=ax.get_zorder()-1, color="white", axis="both")
    ax2.plot(xvals, ydict["CO/O2"], marker="o", zorder=10, color="red")
    ax2.tick_params(axis="y", labelcolor="red")
    ax2.set_ylabel("CO/O"*L"_2"*" ratio", color="red")
    ytix = collect(0:0.2:1.2)
    if maximum(ydict["CO/O2"]) / minimum(ydict["CO/O2"]) >= 10
        ax2.set_yscale("log")
    end
    for side in ["top", "bottom", "left", "right"]
        ax2.spines[side].set_visible(false)
    end

    # past studies
    pastred = "red"
    if pathkey=="water"
        ax2.errorbar(paststudies[lookupkey]["nair"][1], 0.14, capsize=2, zorder=10,
                     yerr=reshape([0.05; 0.43], 2, 1), color=pastred,
                     ecolor=pastred, marker="v")  # Nair 1994 low water
        ax2.errorbar(paststudies[lookupkey]["nair"][2], 0.1, capsize=2, zorder=10,
                     yerr=reshape([0.07; 0.4], 2, 1), color=pastred,
                     ecolor=pastred, marker="^")  # Nair 1994 high water
    else
        ax2.errorbar(paststudies[lookupkey]["nair"], 0.1, capsize=2, zorder=10,
                     yerr=reshape([0.07; 0.4], 2, 1), color=pastred,
                     ecolor=pastred, marker="^")  # Nair 1994 high water
    end
    ax2.scatter(paststudies[lookupkey]["yung"], 0.53, marker="*", s=100,
                color=pastred, zorder=10)  # Yung 1988

    # more ticks for exobase surface temp
    if lookupkey == "temp_exobase"
        ax.set_xticks(ticks=vcat(xvals, [325, 350]))
        [l.set_visible(false) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 2 != 0]
    end

    plotpath = lead*pathkey*"_plots/"*s
    savepath = plotpath*lookupkey*"_CO_and_O2.png"
    savefig(savepath, bbox_inches="tight")
    close(fig)
end

# For the main case in temperatures since they have weird numbers -------------

function makenomdict(abs_or_mr)
    temps = [meanTs, meanTt, meanTe]

    # get the current array
    tfile = "/home/emc/GDrive-CU/Research/Results/TradeoffPlots/Tradeoffs - solar mean/temp_$(meanTsint)_$(meanTtint)_$(meanTeint)/converged_temp_$(meanTsint)_$(meanTtint)_$(meanTeint).h5"
    ncur =  get_ncurrent(tfile)
    N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)
    Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 250e5)
    LA = collect(0e5:2e5:78e5)

    nomTdict = Dict("O2"=>[], "HD"=>[], "HDtop"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                     "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                     "DH"=>[])
    DHprofs = Array{Any}(undef, 1, length(alt)-2)
    # Calculate the things we care about
    # O2 Mixing ratio at surface
    append!(nomTdict["O2"], ncur[:O2][1]/N0)
    # HD mixing ratio
    append!(nomTdict["HD"], ncur[:HD][1]/N0)
    append!(nomTdict["HDtop"], ncur[:HD][end]/Ntop)
    # H2 mixing ratio
    append!(nomTdict["H2"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h) for h in LA]))
    # append!(nomTdict["H2"], sum(ncur[:H2][H2min:H2max]))
    # D Mixing ratio
    append!(nomTdict["D"], ncur[:D][end]/Ntop)
    # H mixing ratio
    append!(nomTdict["H"], ncur[:H][end]/Ntop)
    # D/H at 150 km
    append!(nomTdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
    # CO mixing ratio
    append!(nomTdict["CO"], ncur[:CO][1]/N0)
    # CO/O2 ratio
    append!(nomTdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
    # O3 mixing ratio
    append!(nomTdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
    # H and D fluxes
    Hf = get_H_fluxes(tfile, 1.2e8, temps)
    Df = get_D_fluxes(tfile, 1.2e8, temps)
    append!(nomTdict["Hflux"], Hf)
    append!(nomTdict["Dflux"], Df)
    # fractionation factor 
    append!(nomTdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
    # D/H profile
    DHprofs[1, :] = ncur[:D] ./ ncur[:H]  # altitude profile

    return nomTdict
end

# Main plotting functions ------------------------------------------------------
function make_water_plots(water_x, d, DHdata, q, nom_i, s)
    #=
    Makes the plots for water variation experiment
    water_x: a list of the water mixing ratios, in string format
    d: a dictionary of values of different measurables as function of
       water mixing ratio in lower atmosphere
    DH data: vertical profiles of D/H by experiment
    q: " absolute" or " mixing ratio", just used for labeling
    nom_i: index in water_x of nominal value of water mixing ratio
    s: either "abs" or "mr" 
    =#
    normed_dict = Dict()

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = "Water"*q*" (ppm)"

    # Do the individual plots
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HDtop", "OD", "DO2", "HDO2"], [k for k in keys(d)])
    for i in plots
        # set up plot
        fig, ax = subplots(figsize=(6,4))
        ax.set_facecolor("#ededed")
        grid(zorder=0, color="white")
        xscale("log")
        for side in ["top", "bottom", "left", "right"]
            ax.spines[side].set_visible(false)
        end
        plot(water_x, d[i], marker="o", zorder=10) 
        ax.axvline(nom_i, color="black", label="Nominal value")

        # past studies
        if i=="f"
            scatter(15, 0.32, marker="*", s=100, color="xkcd:tangerine", zorder=10)  # Yung 1988
        end
        
        # set axes labels
        xlabel(L"Total atmospheric water content (pr $\mu$m)")
        ylabdict = Dict("DH"=>"D/H ratio (/1.6e-4) @ 150 km",
                        "D"=>"D"*q*" at exobase",
                        "H"=>"H"*q*" at exobase",
                        "Dflux"=>L"$\phi_D$ (cm$^{-2}$s$^{-1}$)",
                        "Hflux"=>L"$\phi_H$ (cm$^{-2}$s$^{-1}$)",
                        "f"=>"fractionation factor (f)",
                        "O2"=>"O"*L"_2"*q,
                        "HD"=>"HD"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>L"O$_3$ (#/cm$^{-2}$)")
        ylabel(ylabdict[i])

        # Certain measureables have a range not suited to log scale
        nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
        if ~(i in nologplease)#maximum(odict[i])/minimum(odict[i]) > 10#
            yscale("log")
            if i in ["Hflux", "HD", "H", "H2", "O3", "O2"]
                ylim(minimum(d[i])/2, maximum(d[i])*2)
            end
        end

        # set savepath
        plotpath = lead*"water_plots/" * s
        savepath = "water_"*i*".png"

        savefig(plotpath*savepath, bbox_inches="tight")
        close(fig)
    end

    # make special plots
    CO_O2_plot(water_x, d, xlab, "water", nom_i, s)
    DH_alt_prof_plot(DHdata, water_x, "water", s)
end

function make_Oflux_plots(phiO, phiO_str, d, DHdata, q, nom_i, s)
    #=
    Makes the plots for O flux variation experiment
    phiO: a list of the O flux values
    phiO_str: same thing but strings, used for plotting
    d: a dictionary of values of different measurables as function of
       O flux at exobase
    DHdata: altitude profiles of D/H by experiment.
    q: " absolute" or " mixing ratio", just used for labeling
    nom_i: index in phiO of nominal value of O flux
    s: either "abs" or "mr" 
    =#
    normed_dict = Dict()

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = L"$\phi_O$ (cm$^{-2}$s$^{-1}$)"

    # Make individual plots showing change 
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HDtop", "OD", "DO2", "HDO2"], [k for k in keys(d)])
    for i in plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        ax.set_facecolor("#ededed")
        grid(zorder=0, color="white")
        for side in ["top", "bottom", "left", "right"]
            ax.spines[side].set_visible(false)
        end
        plot(phiO_str, d[i], marker="o", zorder=10) 
        ax.axvline(nom_i, color="black", label="Nominal value")

        # past studies
        if i=="f"
            scatter(findfirst(isequal(8e7), phiO)-1, 0.32, marker="*", s=100,
                    color="green", zorder=10)  # Yung 1988
        end
        
        # set axes labels
        xlabel(L"$\phi_O$ (cm$^{-2}$s$^{-1}$)",)
        xticks(rotation=45)
        ylabdict = Dict("DH"=>"D/H ratio (/1.6e-4) @ 150 km",
                        "D"=>"D"*q*"at exobase",
                        "H"=>"H"*q*"at exobase",
                        "Dflux"=>L"$\phi_D$ (cm$^{-2}$s$^{-1}$)",
                        "Hflux"=>L"$\phi_H$ (cm$^{-2}$s$^{-1}$)",
                        "f"=>"fractionation factor (f)",
                        "O2"=>"O"*L"_2"*q,
                        "HD"=>"HD"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>"O"*L"_3"*q)
        ylabel(ylabdict[i])

        # Certain measureables have a range not suited to log scale
        nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
        if ~(i in nologplease)#maximum(d[i])/minimum(d[i]) > 10#
            yscale("log")
            if i in ["Hflux", "HD", "O3"]
                ylim(minimum(d[i])/2, maximum(d[i])*2)
            end
        end

        # set savepath
        plotpath = lead*"Oflux_plots/"*s
        savepath = plotpath*"O_flux_"*i*".png"
        
        savefig(savepath, bbox_inches="tight")
        close(fig)
    end
    
    # make special plots
    CO_O2_plot(phiO_str, d, xlab, "Oflux", nom_i, s)
    DH_alt_prof_plot(DHdata, phiO_str, "Oflux", s)
end

function make_T_plots(T, T_str, d, DHdata, exp, q, nomT, s)
    #=
    Makes the plots for temperature variation experiment
    d: a dictionary of values of different measurables as function of
       temperatures at 3 points in atmosphere: surface, tropopause, exobase
    DHdata: altitude profiles of D/H by experiment.
    T: a list of the temperature values on the x axis
    T_str: same thing but strings, used for plotting
    q: " absolute" or " mixing ratio", just used for labeling
    nomT: a dictionary of nominal values of temperature at the 3 points
    s: either "abs" or "mr" 
    =#

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22
    
    xlab = exp*" temperature (K)"

    # loop through the parameters of interest and plot them
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HDtop", "OD", "DO2", "HDO2"], [k for k in keys(d)])
    for i in plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        better_plot_bg(ax)
        plot(T[exp], d[i], marker="o", zorder=10) 
        ax.axvline(nomT[exp], color="black", label="Nominal value")
        xlabel(xlab)
        xticks(ticks=T[exp], labels=T_str[exp], rotation=45)
        ylabdict = Dict("DH"=>"D/H ratio (/1.6e-4) @ 150 km",
                        "D"=>"D"*q*"at exobase",
                        "H"=>"H"*q*"at exobase",
                        "Dflux"=>L"$\phi_D$ (cm$^{-2}$s$^{-1}$)",
                        "Hflux"=>L"$\phi_H$ (cm$^{-2}$s$^{-1}$)",
                        "f"=>"fractionation factor (f)",
                        "O2"=>"O"*L"_2"*q,
                        "HD"=>"HD"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>"O"*L"_3"*q)
        ylabel(ylabdict[i])

        # past studies
        if exp=="surface"
            if i == "f"
                scatter(220, 0.32, marker="d", color="xkcd:tangerine", 
                        zorder=10)  # Yung 1988
            end
        elseif exp=="tropopause"
            if i == "f"
                scatter(125, 0.055, marker="v", color="xkcd:golden", 
                        zorder=10)  #Kras 2002 solar minimum
                scatter(125, 0.082, marker="s", color="xkcd:blood orange", 
                        zorder=10)  #Kras 2002 solar mean
                scatter(125, 0.167, marker="^", color="xkcd:scarlet", 
                        zorder=10)  #Kras 2002 solar max
                errorbar(140, 0.14, yerr=reshape([0.05; 0.43], 2, 1), 
                         fmt="s", color="purple", ecolor="gray",
                         zorder=10, capsize=2)  # Nair 1994 low water
            end
        elseif exp=="exobase"
            if i == "f"
                scatter(200, 0.055, marker="v", color="xkcd:golden", 
                        zorder=10)  #Kras 2002 solar minimum
                scatter(270, 0.082, marker="s", color="xkcd:blood orange", 
                        zorder=10)  #Kras 2002 solar mean
                scatter(350, 0.167, marker="^", color="xkcd:scarlet", 
                        zorder=10)  #Kras 2002 solar max
                errorbar(288, 0.14, yerr=reshape([0.05; 0.43], 2, 1), 
                         fmt="s", color="xkcd:hunter green", 
                         ecolor="xkcd:dark sage", zorder=10, capsize=2)  # Nair 1994 low water
                xticks(ticks=vcat(T[exp], [325, 350]), 
                       labels=vcat(T_str[exp], 
                       ["325", "350"]), rotation=45)
            end
        end 

        # Certain measureables have a range not suited to log scale
        nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
        if ~(i in nologplease)#maximum(d[i])/minimum(d[i]) > 10#
            yscale("log")
            if i in ["Hflux", "HD"]
                ylim(minimum(d[i])/2, maximum(d[i])*2)
            end
        end

        # set savepath
        plotpath = lead*"temp_plots/"*s
        savepath = plotpath*join(["temp", exp], "_")*"_"*i*".png"
        savefig(savepath, bbox_inches="tight")
        close(fig)
    end

    # make special plots
    CO_O2_plot(T[exp], d, exp*" temperature", "temp", nomT[exp], s, "_"*exp)
    DH_alt_prof_plot(DHdata, T_str[exp], "temp", s, "_"*exp, 
                     latexstring("T_{$(exp[1:3])}"))
end

function make_rel_change_plots(output_MR, output_abs, ex, SVP)
    #=
    output_MR: 
    output_abs:
    ex: experiment type, the usual key of "water", "O flux", "surface" etc.
    SVP: whether to plot for experiments where SVP was held constant 
         (SVP="const") to one temp profile, or for experiments where SVP was 
         allowed to follow the temp profile (SVP="vary")
    =#
    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    if SVP=="const"
        surfvals = [150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 
                    230.0, 240.0, 250.0, 260.0, 270.0]
    elseif SVP=="vary"
        surfvals = [190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0]
    end
    
    # the official values of each experiment that I ran 
    # water vals in MR: ["5e-6", "5e-5", "8e-4", "2.2e-3", "4.4e-3", "9e-3", "1.35e-2"]
    xvals = Dict("water"=>[0.1, 1, 10, 25, 50, 100, 150],  # in pr μm, total column
                 "O flux"=>["3e7", "4e7", "5e7", "6e7", "7e7", "8e7", "9e7", 
                            "1.0e8", "1.1e8", "1.2e8", "1.3e8", "1.4e8", "1.5e8", 
                            "1.6e8"],
                 "surface"=>surfvals,
                 "tropopause"=>[70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 
                                140.0, 150.0, 160.0],
                 "exobase"=>[150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 
                             300.0, 325.0, 350.0])

    # something to add onto the x label for each plot
    xlab_addn = Dict("water"=>" mixing ratio (ppm)", "O flux"=>L" (cm$^{-2}s^{-1}$)", 
                     "surface"=>" temperature", "tropopause"=>" temperature", 
                     "exobase"=>" temperature")

    # first plot - compare with data ---------------------------------------------
    fig, ax1 = subplots(figsize=(7,6))
  
    # Data from MAVEN or other missions and styling to apply on all axes
    shade = "gainsboro"
    better_plot_bg(ax1)
    # specific formatting by experiment type
    if ex=="exobase"
        ax1.set_xticks(150:25:350)
        ax1.axvspan(150, 300, color=shade)
    elseif ex=="tropopause"
        ax1.axvspan(75, 160, color=shade)
    elseif ex=="surface"
        ax1.axvspan(180, 270, color=shade)
    elseif ex=="O flux"
        ax1.axvspan(5, 6, color=shade) # doubled to 8.6e7 from 4.3e7 per Mike
        ax1.tick_params(rotation=45, axis="x")
    elseif ex=="water"
        # ax1.axvspan(80, 200, color="xkcd:brick red", alpha=0.3) #dust storm
        ax1.axvspan(8, 12, color=shade) #nominal
        # ax1.axvspan(30, 80, color="#67a9cf", alpha=0.3) # summer pole
    end

    # nominal value plot location or index, 1-indexed for Julia
    nom_i_julia = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs, 
                       "O flux"=>findfirst(isequal("1.2e8"), xvals["O flux"]), 
                       "water"=>10)
    # Passing the values to PyPlot requires 0-indexing so here it is:
    nom_i_py = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs, 
                    "O flux"=>findfirst(isequal("1.2e8"), xvals["O flux"])-1, 
                    "water"=>10)
    
    if ex in ["exobase", "tropopause", "surface"]
        unit = " (K)"
    else
        unit = ""
    end

    ax1.set_xlabel(ex*xlab_addn[ex]*unit)
    ax1.set_ylabel(L"(Model-Obs)/$\sigma$")
    c1 = ["#10007A", "#2F7965", "#e16262", "#D59A07"]  # just colors

    # set the data for comparison. CO MR (Nair+1994), O2 MR (Nair+1994), 
    # H2 abundance (Kras&Feldman 2001), O3 μm-atm (Clancy 2016)
    data = [6e-4, 1.2e-3, 15, 1.18]  # CO MR, O2 MR, H2 ppm, O3 μm-atm
    s = [1.5e-4, 0.2e-3, 5, 0.7]    # sigmas (uncertainty) on each
    
    # here, set a dummy var to handle fact that temperature has one more index 
    # for lookup. 
    if ex in ["water", "O flux"]
        MRdictvar = output_MR
        ABSdictvar = output_abs
    else
        MRdictvar = output_MR[ex]
        ABSdictvar = output_abs[ex]
    end

    # calculate the relativeness of each point wrt data
    COdiff = (MRdictvar["CO"] .- data[1])/s[1]
    O2diff = (MRdictvar["O2"] .- data[2])/s[2]
    # H2diff = (ABSdictvar["H2"] .- data[3])/s[3]  # absolute abundance
    H2diff = (MRdictvar["H2"]./1e-6 .- data[3])/s[3] # ppm
    O3diff = (areadensity_to_micron_atm(ABSdictvar["O3"]) .- data[4])/s[4] 


    # plot the actual stuff   
    sz = 10
    ax_12 = ax1.twinx()
    ax1.plot(xvals[ex], COdiff, marker="o", ms=sz, color=c1[1], 
               zorder=10, label="CO")
    ax1.plot(xvals[ex], O2diff, marker="x", ms=sz, color=c1[2], 
               zorder=10, label=L"O$_2$")
    ax1.plot(xvals[ex], H2diff, marker="*", ms=sz, color=c1[3], 
               zorder=10, label=L"H$_2$")
    ax1.plot(xvals[ex], O3diff, marker="^", ms=sz, color=c1[4], 
               zorder=10, label=L"O$_3$")

    # second axis to show CO/O2
    ax_12.plot(xvals[ex], ABSdictvar["CO"] ./ ABSdictvar["O2"], color="cornflowerblue", marker="D")
    ax_12.plot(xvals[ex], 0.6*fill!(similar(xvals[ex]), 1), color="cornflowerblue", alpha=0.5)
    ax_12.text(xvals[ex][end], 0.6, "0.6", ha="right", va="bottom", color="cornflowerblue")
    ax_12.set_ylabel(L"CO/O$_2$", color="cornflowerblue")
    ax_12.tick_params(axis="y", labelcolor="cornflowerblue", which="both")
    for side in ["top", "bottom", "left", "right"]
        ax_12.spines[side].set_visible(false)
    end


    # plot the nominal values
    medgray = "#444444"
    ax1.axvline(nom_i_py[ex], color=medgray, zorder=5)
    ax1.axhline(0, color="black")

    # Text on plots 
    if ex=="exobase"
        ax1.text(208, 5, "Nominal\ncase", color=medgray, ha="left")
        ax1.text(150, -2.6, "CO", color=c1[1], ha="left", va="top")
        ax1.text(170, -4.8, L"O$_2$", color=c1[2], ha="left", va="top")
        ax1.text(160, 7, L"H$_2$", color=c1[3], ha="left", va="top")
        ax1.text(150, -0.1, L"O$_3$", color=c1[4], ha="left", va="top")
    elseif ex=="O flux"
        ax1.text(8.7, 1, "Nominal\ncase", color=medgray, ha="right")
    elseif ex=="surface"
        ax1.text(218, 30, "Nominal\ncase", color=medgray, ha="left")
        ax1.text(150, -5.2, "CO", color=c1[1], ha="left", va="top")
        ax1.text(145, 25, L"O$_2$", color=c1[2], ha="left", va="top")
        ax1.text(150, 5, L"H$_2$", color=c1[3], ha="left", va="top")
        ax1.text(155, 45, L"O$_3$", color=c1[4], ha="left", va="top")
        ax1.set_yticks(-10:10:50)
        ax1.set_xticks(150:20:280)
        ax1.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
        ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
        ax1.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    elseif ex=="tropopause"
        ax1.text(109, 2, "Nominal\ncase", color=medgray, ha="left")
        ax1.text(70, -3, "CO", color=c1[1], ha="left", va="top")
        ax1.text(90, -4.2, L"O$_2$", color=c1[2], ha="left", va="top")
        ax1.text(70, 3.2, L"H$_2$", color=c1[3], ha="left", va="top")
        ax1.text(70, -0.5, L"O$_3$", color=c1[4], ha="left", va="top")
        ax1.set_xticks(70:20:160)
        ax1.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
        ax1.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
    elseif ex=="water"
        ax1.set_xscale("log")
        # ax1.text(81, 4, " Dust \nstorm", color="xkcd:brick red")
        # ax1.text(25, 6.5, "Summer\n   pole", color="navy")
        ax1.text(9, 10, "Nominal\ncase", color=medgray, ha="right")
        ax1.text(0.2, 11, "CO", color=c1[1], ha="left", va="top")
        ax1.text(2.5, -1, L"O$_2$", color=c1[2], ha="left", va="top")
        ax1.text(4, 2, L"H$_2$", color=c1[3], ha="left", va="top")
        ax1.text(0.1, 3, L"O$_3$", color=c1[4], ha="left", va="top")
    end
    #ax1.legend(bbox_to_anchor=(1.05, 1))
    savefig(lead*"metrics_tradeoff_"*ex*".png", bbox_inches="tight")


    # second plot: H2, HD, H, D, Hflux, Dflux ----------------------------------
    fig, ax2 = subplots(figsize=(7,6))
   
    # Data from MAVEN or other missions and styling to apply on all axes
    shade = "gainsboro"
    better_plot_bg(ax2)
    # specific formatting by experiment type
    if ex=="exobase"
        ax2.set_xticks(150:25:350)
        ax2.axvspan(150, 300, color=shade)
    elseif ex=="tropopause"
        ax2.axvspan(75, 160, color=shade)
    elseif ex=="surface"
        ax2.axvspan(180, 270, color=shade)
    elseif ex=="O flux"
        ax2.axvspan(5, 6, color=shade) # doubled to 8.6e7 from 4.3e7 per Mike
        ax2.tick_params(rotation=45, axis="x")
    elseif ex=="water"
        # ax2.axvspan(80, 200, color="xkcd:brick red", alpha=0.3) #dust storm
        ax2.axvspan(8, 12, color=shade) #nominal
        # ax2.axvspan(30, 80, color="#67a9cf", alpha=0.3) # summer pole
    end

    ax2.set_ylabel(L"X/X$_{nominal}$")
    c2 = ["#6270d9","#95b034","#d14f58","#5ca85c","#ce6d30"]

    # to do the division/normalization we need to re-find the index of the value
    # against which to normalize. To do this, look in the xvals by experiment, 
    # then index that according to whether we cut it or not. There doesn't seem 
    # to be a cleaner way to do the indexing of cut.
    normidx = Dict("exobase"=>findfirst(isequal(meanTe), xvals["exobase"][1:1:length(xvals["exobase"])]), 
                   "tropopause"=>findfirst(isequal(meanTt), xvals["tropopause"][1:1:length(xvals["tropopause"])]), 
                   "surface"=>findfirst(isequal(meanTs), xvals["surface"][1:1:length(xvals["surface"])]),
                   "water"=>findfirst(isequal(10), xvals["water"][1:1:length(xvals["water"])]),
                   "O flux"=>findfirst(isequal("1.2e8"), xvals["O flux"][1:1:length(xvals["O flux"])]))

    if ex in ["water", "O flux"]
        HDdiff = normalize(MRdictvar["HD"], normidx[ex])
        Hdiff = normalize(MRdictvar["H"], normidx[ex])
        Ddiff = normalize(MRdictvar["D"], normidx[ex])
        Hfdiff = normalize(MRdictvar["Hflux"], normidx[ex])
        Dfdiff = normalize(MRdictvar["Dflux"], normidx[ex])
    else
        nomdict = makenomdict("mr") # get info for nom case which has values not % by 10
        HDdiff = normalize_val(MRdictvar["HD"], nomdict["HD"])
        Hdiff = normalize_val(MRdictvar["H"], nomdict["H"])
        Ddiff = normalize_val(MRdictvar["D"], nomdict["D"])
        Hfdiff = normalize_val(MRdictvar["Hflux"], nomdict["Hflux"])
        Dfdiff = normalize_val(MRdictvar["Dflux"], nomdict["Dflux"])
    end

    # plot the actual stuff
    ax2.plot(xvals[ex], HDdiff, marker="o", ms=sz, color=c2[1], zorder=10, label="HD")
    ax2.plot(xvals[ex], Hdiff, marker="x", ms=sz, color=c2[2], zorder=10, label="H")
    ax2.plot(xvals[ex], Ddiff, marker="*", ms=sz, color=c2[3], zorder=10, label="D")
    ax2.plot(xvals[ex], Hfdiff, marker="^", ms=sz, color=c2[4], zorder=10, label=L"\phi_H")
    ax2.plot(xvals[ex], Dfdiff, marker="D", ms=sz, color=c2[5], zorder=10, label=L"\phi_D")

    # Text on plots 
    if ex=="exobase"
        ax2.text(208, 100, "Nominal\ncase", color=medgray, ha="left", va="top")
        ax2.text(150, 5, "HD", color=c2[1], ha="left", va="top")
        ax2.text(175, 3, "H", color=c2[2], ha="left", va="top")
        ax2.text(150, 0.4, "D", color=c2[3], ha="left", va="top")
        ax2.text(325, 2, L"$\phi_H$", color=c2[4], ha="left", va="top")
        ax2.text(150, 0.15, L"$\phi_D$", color=c2[5], ha="left", va="top")
    # elseif ex=="O flux"
    #     ax2.text(8.7, 1, "Nominal\ncase", color=medgray, ha="right")
    elseif ex=="surface"
        ax2.text(219, 1.6, "Nominal\ncase", color=medgray, ha="left", va="top")
        ax2.text(155, 0.6, "HD", color=c2[1], ha="left", va="top")
        ax2.text(147, 1.35, "H", color=c2[2], ha="left", va="top")
        ax2.text(155, 1.6, "D", color=c2[3], ha="left", va="top")
        ax2.text(255, 0.95, L"$\phi_H$", color=c2[4], ha="left", va="top")
        ax2.text(255, 1.2, L"$\phi_D$", color=c2[5], ha="left", va="top")
        ax2.set_xticks(150:20:280)
        ax2.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
        # ax2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
        # ax2.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    elseif ex=="tropopause"
        ax2.text(110, 3.25, "Nominal\ncase", color=medgray, ha="left", va="top")
        ax2.text(70, 2.05, "HD", color=c2[1], ha="left", va="top")
        ax2.text(70, 1.3, "H", color=c2[2], ha="left", va="top")
        ax2.text(70, 0.8, "D", color=c2[3], ha="left", va="top")
        ax2.text(150, 1.25, L"$\phi_H$", color=c2[4], ha="left", va="top")
        ax2.text(70, 0.4, L"$\phi_D$", color=c2[5], ha="left", va="top")
        ax1.set_xticks(70:20:160)
        ax1.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
        ax1.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
    elseif ex=="water"
        ax2.set_xscale("log")
        ax2.text(11, 0.9, "Nominal\ncase", color=medgray, ha="left", va="top")
        ax2.text(30, 0.965, "HD", color=c2[1], ha="left", va="top")
        ax2.text(0.1, 0.982, "H", color=c2[2], ha="left", va="top")
        ax2.text(0.1, 0.91, "D", color=c2[3], ha="left", va="top")
        ax2.text(30, 0.995, L"$\phi_H$", color=c2[4], ha="left", va="top")
        ax2.text(0.12, 0.95, L"$\phi_D$", color=c2[5], ha="left", va="top")
    end


    # Other species that can be plotted, or species at other locations
    # HDtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ODtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # HDO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # DO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ax2.plot(xvals[ex], HDtopdiff, marker="o", color="#7D3403", zorder=10)
    # ax2.plot(xvals[ex], ODtopdiff, marker="o", color="red", zorder=10)
    # ax2.plot(xvals[ex], HDO2topdiff, marker="o", color="blue", zorder=10)
    # ax2.plot(xvals[ex], DO2topdiff, marker="o", color="purple", zorder=10)

    # nominal value
    ax2.axvline(nom_i_py[ex], color=medgray, zorder=5)
    ax2.axhline(1, color=medgray, zorder=5)
    
    if ex in ["exobase", "tropopause", "surface"]
        unit = " (K)"
        if ex=="exobase"
            ax2.set_yscale("log")
        end
    else
        unit = ""
    end
    ax2.set_xlabel(ex*xlab_addn[ex]*unit)
    #ax2.legend(bbox_to_anchor=(1.05, 1))
    savefig(lead*"compare_nominal_"*ex*".png", bbox_inches="tight")


    # third plot: f ------------------------------------------------------------
    fig, ax = subplots(figsize=(7,6))
  
    # Data from MAVEN or other missions and styling to apply on all axes
    shade = "gainsboro"
    better_plot_bg(ax)
    # specific formatting by experiment type
    if ex=="exobase"
        ax.set_xticks(150:25:350)
        ax.axvspan(150, 300, color=shade)
        ax.text(203, 0.01, "Nominal\ncase", color=medgray, ha="right")
    elseif ex=="tropopause"
        ax.axvspan(75, 160, color=shade)
        ax.text(110, 0.005, "Nominal\ncase", color=medgray, ha="left")
    elseif ex=="surface"
        ax.set_xticks(150:20:270)
        ax.axvspan(180, 270, color=shade)
        ax.text(218, 2.5e-3 , "Nominal\ncase", color=medgray)
        ax.grid(zorder=0, color="white", which="both")
    elseif ex=="O flux"
        ax.axvspan(5, 6, color=shade) # doubled to 8.6e7 from 4.3e7 per Mike
        ax.tick_params(rotation=45, axis="x")
        ax.text(8.7, 0.00112, "Nominal\ncase", color=medgray, ha="right")
    elseif ex=="water"
        # ax.axvspan(80, 200, color="xkcd:brick red", alpha=0.3) #dust storm
        ax.axvspan(8, 12, color=shade) #nominal
        # ax.axvspan(30, 80, color="#67a9cf", alpha=0.3) # summer pole
        ax.text(9, 0.00095, "Nominal\ncase", color=medgray, ha="right", va="top")
        # ax.text(50, 0.001025, "Summer\npole", color="navy", ha="center", va="top")
        # ax.text(150, 0.001, "Dust\nstorm", color="xkcd:brick red", ha="center", va="top")
        ax.set_xscale("log")
    end

    ax.set_ylabel(L"$f$", color="black")
    ax.set_xlabel(ex*xlab_addn[ex])
    if ex in ["water", "Oflux"]
        linecol = "cornflowerblue"
    else
        linecol = "black"
    end
    ax.plot(xvals[ex], MRdictvar["f"], marker="o", color=linecol, zorder=10)
    ax.axvline(nom_i_py[ex], color=medgray, zorder=5)

    ylims = Dict("tropopause"=>[1e-4, 1e-2],
                 "exobase"=>[1e-5, 2e-1],
                 "surface"=>[1e-3, 3e-3])

    if ex in ["exobase", "tropopause", "surface"]  # in these cases, we'll have 2 axes.
        # some adjustments to the main axes, which shows raw value
        ax.set_yscale("log") 
        ax.set_ylim(ylims[ex][1], ylims[ex][2])
        ax.set_xlabel(ex*xlab_addn[ex]*" (K)")

        # plot the secondary axis showing increase
        fdiff = normalize_val(MRdictvar["f"], nomdict["f"])
        ax_2 = ax.twinx()
        ax_2.plot(xvals[ex], fdiff, marker="o", color="xkcd:vivid purple")
        ax_2.set_ylabel(L"$f$/$f_{nominal}$", color="xkcd:vivid purple")
        ax_2.tick_params(axis="y", labelcolor="xkcd:vivid purple", which="both")
        for side in ["top", "bottom", "left", "right"]
            ax_2.spines[side].set_visible(false)
        end
    end
    
    savefig(lead*"f_tradeoff_"*ex*".png", bbox_inches="tight")
end

# Analyzation functions (main routines) ----------------------------------------

function analyze_water(abs_or_mr, allDbearers, make_plots=false, path=lead)
    # Establish parameters, filenames, etc
    watervals = [0.1, 1, 10, 25, 50, 100, 150]
    watervals_str = ["1.48e-6", "2.1e-5", "8.1e-4", "2.18e-3", "4.46e-3", "9.02e-3", "1.358e-2"]
    wfilelist = [path*"water_"*w*"/converged_water_"*w*".h5" for w in watervals_str]
    temps = [meanTs, meanTt, meanTe]
    oflux = 1.2e8
    q = abs_or_mr == "abs" ? " abundance" : " mixing ratio" # for labels
    mean_idx = findfirst(isequal(10), watervals) - 1 
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"
    
    # Establish variables to store data on simulations
    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(watervals), length(alt)-2)  
    wdict = Dict{String, Array}("O2"=>[], "HD"=>[], "HDtop"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                 "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                 "DH"=>[])
    if allDbearers
        wdict["OD"]=>[]
        wdict["HDO2"]=>[]
        wdict["DO2"]=>[]
    end

    # loop water files and collect data 
    i = 1
    for wfile in wfilelist
        # get the current array
        ncur = get_ncurrent(wfile)

        N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)
        Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 250e5)
        LA = collect(0e5:2e5:78e5)
        # H2min = n_alt_index[140e5] # Kras & Feldman 2001 measured H2 above 140 km
        # H2max = n_alt_index[200e5]

        # Calculate the things we care about
        # O2 Mixing ratio at surface
        append!(wdict["O2"], ncur[:O2][1]/N0)
        # HD mixing ratio
        append!(wdict["HD"], ncur[:HD][1]/N0)
        append!(wdict["HDtop"], ncur[:HD][end]/Ntop)
        # H2 mixing ratio
        append!(wdict["H2"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h) for h in LA])) # ppm in lower atmo
        # append!(wdict["H2"], sum(ncur[:H2][H2min:H2max]))
        # D Mixing ratio
        append!(wdict["D"], ncur[:D][end]/Ntop)
        # H mixing ratio
        append!(wdict["H"], ncur[:H][end]/Ntop)
        # D/H at 150 km
        append!(wdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # CO mixing ratio
        append!(wdict["CO"], ncur[:CO][1]/N0)
        # CO/O2 ratio
        append!(wdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 mixing ratio
        append!(wdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
        # H and D fluxes
        Hf = get_H_fluxes(wfile, oflux, temps)
        Df = get_D_fluxes(wfile, oflux, temps)
        append!(wdict["Hflux"], Hf)
        append!(wdict["Dflux"], Df)
        # fractionation factor 
        append!(wdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
        # Other D bearing species???
        if allDbearers
            append!(wdict["OD"], ncur[:OD][end]/Ntop)
            append!(wdict["HDO2"], ncur[:HDO2][end]/Ntop)
            append!(wdict["DO2"], ncur[:DO2][end]/Ntop)
        end
        i += 1
    end

    if make_plots == true
        make_water_plots(watervals, wdict, DHprofs, q, mean_idx, subfolder)
    end

    println("Finished water plots")
    return wdict
end

function analyze_Oflux(abs_or_mr, allDbearers, make_plots=false, path=lead)
    # Establish important parameters, files, etc
    Ofluxvals = [3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7, 1e8, 1.1e8, 1.2e8, 1.3e8, 
                 1.4e8, 1.5e8, 1.6e8]
    Ofluxvals_str = ["3e7", "4e7", "5e7", "6e7", "7e7", "8e7", "9e7", "1.0e8", 
                     "1.1e8", "1.2e8", "1.3e8", "1.4e8", "1.5e8", "1.6e8"]
    Ofilelist = [path*"Oflux_"*o*"/converged_Oflux_"*o*".h5" 
                 for o in Ofluxvals_str]
    temps = [meanTs, meanTt, meanTe]
    mean_idx = findfirst(isequal(1.2e8), Ofluxvals) - 1
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # Establish variables to store data on simulations
    odict = Dict("O2"=>[], "HD"=>[], "HDtop"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                 "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                 "DH"=>[])
    if allDbearers
        odict["OD"]=>[]
        odict["HDO2"]=>[]
        odict["DO2"]=>[]
    end

    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(Ofluxvals_str), length(alt)-2)

    # loop through O flux files 
    i = 1
    for (oflux, ofile) in zip(Ofluxvals, Ofilelist)
        # calculate the things we care about
        # get the current array
        ncur = get_ncurrent(ofile)

        N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)
        Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 250e5)
        LA = collect(0e5:2e5:78e5)
        # H2min = n_alt_index[140e5]
        # H2max = n_alt_index[200e5]
        
        # Calculate the things we care about
        # O2 Mixing ratio at surface
        append!(odict["O2"], ncur[:O2][1]/N0)
        # HD mixing ratio
        append!(odict["HD"], ncur[:HD][1]/N0)
        append!(odict["HDtop"], ncur[:HD][end]/Ntop)
        # H2 mixing ratio
        append!(odict["H2"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h) for h in LA]))
        # append!(odict["H2"], sum(ncur[:H2][H2min:H2max]))
        # D Mixing ratio
        append!(odict["D"], ncur[:D][end]/Ntop)
        # H mixing ratio
        append!(odict["H"], ncur[:H][end]/Ntop)
        # D/H at 150 km
        append!(odict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # CO mixing ratio
        append!(odict["CO"], ncur[:CO][1]/N0)
        # CO/O2 ratio
        append!(odict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 mixing ratio
        append!(odict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
        # H and D fluxes
        Hf = get_H_fluxes(ofile, oflux, temps)
        Df = get_D_fluxes(ofile, oflux, temps)
        append!(odict["Hflux"], Hf)
        append!(odict["Dflux"], Df)
        # fractionation factor 
        append!(odict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profiles
        i += 1
        # Other D bearing species???
        if allDbearers
            append!(odict["OD"], ncur[:OD][end]/Ntop)
            append!(odict["HDO2"], ncur[:HDO2][end]/Ntop)
            append!(odict["DO2"], ncur[:DO2][end]/Ntop)
        end
    end

    if make_plots == true
        make_Oflux_plots(Ofluxvals, Ofluxvals_str, odict, DHprofs, q, mean_idx, 
                         subfolder)
    end

    println("Finished O flux")
    return odict
end
    
function analyze_T(abs_or_mr, allDbearers, SVP, make_plots=false, path=lead)
    #=
    abs_or_mr: whether using the absolute abundance file or mixing ratio file
    SVP: whether SVP was held constant ("const") for one temp profile or varied
         with temp profile ("vary")
    make_plots: whether to generate the plots
    path: where to save plots
    =#

    if SVP=="const"
        surfvals = [150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 
                    230.0, 240.0, 250.0, 260.0, 270.0]
    elseif SVP=="vary"
        surfvals = [190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 
                    270.0]
    end

    # Set up parameters, filenames, etc
    tvals = Dict("surface"=>surfvals,
                 "tropopause"=>[70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 
                                140.0, 150.0, 160.0],
                 "exobase"=>[150.0, 175.0, 200.0, 225.0, 250.0, 
                             275.0, 300.0, 325.0, 350.0])
    # Dictionaries to store each experiment's data so it can be returned
    all_tdicts = Dict()

    tvals_str = Dict()
    tempfilelist = Dict()
    for k in keys(tvals)
        tvals_str[k] = [string(trunc(Int, x)) for x in tvals[k]]
        if k == "surface"
            tempfilelist[k] = [path*"temp_"*t*"_$(meanTtint)_$(meanTeint)"*"/converged_temp_"*t*"_$(meanTtint)_$(meanTeint).h5" for t in tvals_str[k]]  # TODO: revert if necessary
        elseif k == "tropopause"
            tempfilelist[k] = [path*"temp_$(meanTsint)_"*t*"_$(meanTeint)"*"/converged_temp_$(meanTsint)_"*t*"_$(meanTeint).h5" for t in tvals_str[k]]
        elseif k == "exobase"
            tempfilelist[k] = [path*"temp_$(meanTsint)_$(meanTtint)_"*t*"/converged_temp_$(meanTsint)_$(meanTtint)_"*t*".h5" for t in tvals_str[k]]
        end
    end
    meanT = Dict("surface"=>meanTs, "tropopause"=>meanTt, "exobase"=>meanTe)  # nominal 
    oflux = 1.2e8
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # loop through which temp is varied and construct a list of datapoints
    for experiment in keys(tvals) # loop across the dictionary
        println("Now doing temperature ", experiment)
        tdict = Dict("O2"=>[], "HD"=>[], "HDtop"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                     "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                     "DH"=>[])
        if allDbearers
            tdict["OD"]=>[]
            tdict["HDO2"]=>[]
            tdict["DO2"]=>[]
        end
        
        DHprofs = Array{Any}(undef, length(tvals_str[experiment]), length(alt)-2)

        # now loop through the values for each varied temp
        i = 1
        for (tv, tfile) in zip(tvals[experiment], tempfilelist[experiment])
            # set the temperature profile
            if experiment == "surface" 
                temps = [tv, meanTt, meanTe] 
            elseif experiment == "tropopause"
                temps = [meanTs, tv, meanTe]
            elseif experiment == "exobase"
                temps = [meanTs, meanTt, tv]
            end

            # get the current array
            ncur = get_ncurrent(tfile)
            N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)
            Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 250e5)
            LA = collect(0e5:2e5:78e5)
            # H2min = n_alt_index[140e5]
            # H2max = n_alt_index[200e5]

            # Calculate the things we care about
            # O2 Mixing ratio at surface
            append!(tdict["O2"], ncur[:O2][1]/N0)
            # HD mixing ratio
            append!(tdict["HD"], ncur[:HD][1]/N0)
            append!(tdict["HDtop"], ncur[:HD][end]/Ntop)
            # H2 mixing ratio
            append!(tdict["H2"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h) for h in LA]))
            # append!(tdict["H2"], sum(ncur[:H2][H2min:H2max]))
            # D Mixing ratio
            append!(tdict["D"], ncur[:D][end]/Ntop)
            # H mixing ratio
            append!(tdict["H"], ncur[:H][end]/Ntop)
            # D/H at 150 km
            append!(tdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
            # CO mixing ratio
            append!(tdict["CO"], ncur[:CO][1]/N0)
            # CO/O2 ratio
            append!(tdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
            # O3 mixing ratio
            append!(tdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
            # H and D fluxes
            Hf = get_H_fluxes(tfile, oflux, temps)
            Df = get_D_fluxes(tfile, oflux, temps)
            append!(tdict["Hflux"], Hf)
            append!(tdict["Dflux"], Df)
            # fractionation factor 
            append!(tdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
            # D/H profile
            DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
            i += 1
            # Other D bearing species???
            if allDbearers
                append!(tdict["OD"], ncur[:OD][end]/Ntop)
                append!(tdict["HDO2"], ncur[:HDO2][end]/Ntop)
                append!(tdict["DO2"], ncur[:DO2][end]/Ntop)
            end
        end

        if make_plots == true
            make_T_plots(tvals, tvals_str, tdict, DHprofs, experiment, q, meanT, 
                         subfolder)
        end
        
        println("Finished temperature ", experiment)

        all_tdicts[experiment] = tdict
    end
    println("Finished temps")
    return all_tdicts
end

# Define mean temperatures =====================================================
global meanTs = 216.0
global meanTt = 108.0
global meanTe = 205.0

global meanTsint = 216
global meanTtint = 108
global meanTeint = 205

# Setup and folder Location ====================================================
S = "const" # TODO: change or check each time 
# lead = "/data/GDrive-CU/Research/Results/TradeoffPlots/Tradeoffs - solar mean/"
lead = "/home/emc/GDrive-CU/Research/Results/TradeoffPlots/Tradeoffs - solar mean/"
makeplots = false   # whether to make the plots that go with each experiment
other_deuterated = false
write_new_files = false  # set to true if running for first time after new simulations

# Function calls ===============================================================

water_data_mr = analyze_water("mr", other_deuterated, makeplots)
water_data_abs = analyze_water("abs", other_deuterated, makeplots)

println()

# o_data_abs = analyze_Oflux("abs", other_deuterated, makeplots)
# o_data_mr = analyze_Oflux("mr", other_deuterated, makeplots)

println()

T_data_mr = analyze_T("mr", other_deuterated, S, makeplots)
T_data_abs = analyze_T("abs", other_deuterated, S, makeplots)

if write_new_files
    wd_mr = jldopen(lead*"water_MR_data.jld", "w")
    @write wd_mr water_data_mr
    close(wd_mr)
    wd_abs = jldopen(lead*"water_abs_data.jld", "w")
    @write wd_abs water_data_abs
    close(wd_abs)

    # O_mr = jldopen(lead*"O_MR_data.jld", "w")
    # @write O_mr o_data_mr
    # close(O_mr)
    # O_abs = jldopen(lead*"O_abs_data.jld", "w")
    # @write O_abs o_data_abs
    # close(O_abs)

    T_mr = jldopen(lead*"T_MR_data.jld", "w")
    @write T_mr T_data_mr
    close(T_mr)
    T_abs = jldopen(lead*"T_abs_data.jld", "w")
    @write T_abs T_data_abs
    close(T_abs)
end

# Relative changes plot with two panels ========================================
# panel 1: metrics CO, O2, O3, and H2
# panel 2: all the other measureablesbles except fractionation factor
# panel 3: frationation factor

make_rel_change_plots(water_data_mr, water_data_abs, "water", S)
# make_rel_change_plots(o_data_mr, o_data_abs, "O flux", S)
make_rel_change_plots(T_data_mr, T_data_abs, "surface", S)
make_rel_change_plots(T_data_mr, T_data_abs, "tropopause", S)
make_rel_change_plots(T_data_mr, T_data_abs, "exobase", S)