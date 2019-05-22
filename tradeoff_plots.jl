################################################################################
# make_tradeoff_plots.jl - make trade off plots showing H and D flux and other
# variables as a function of random things such as water column. 

# Eryn Cangi
# 5 April 2019
# Currently tested for Julia: 0.7
################################################################################
using PyPlot
using HDF5
using LaTeXStrings
using PyCall
using PlotUtils

# Functions ====================================================================
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

function get_H_fluxes(readfile, oflux, temps)
    #=
    Produces an array of H and D fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
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

function get_D_fluxes(readfile, oflux, temps)
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
   =#
   R_HDO = 434.8 # J/kgK
   L_HDO = LHDO(T) * 1000 / 0.019 # kJ/mol * J/kJ * mol / kg = J / kg
   T0 = 373.85 # boiling temp for liquid HDO
   P0 = 101325 # standard pressure
   return P0*exp(-(L_HDO/R_HDO)*(1/T - 1/T0))
end

Psat(T::Float64) = (133.3/(10^6 * boltzmannK * T))*(10^(-2445.5646/T 
                    + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 
                    - 6.757169))

function speciesbcs_VARIABLE(species, oflux, temps)

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

# Special functions for this file ==============================================

searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function normalize(arr, base_i)
    normed_arr = arr ./ arr[base_i]
    return normed_arr
end

function DH_HD_alt_prof_plot(DHproflist, HDproflist, experiments, varyable, optext="", optlegend="")
    #=
    DHproflist: An array with D/H profiles in each row
    HDproflist: same as DHproflist but for the species HD
    experiments: a list of experiment identifiers as strings
    varyable: water, Oflux, or temp (for putting in correct folder)
    optext: an optional extension that will be appended to varyable, for 
            the temperature case where it goes in the temp_tradeoff_plots 
            folder but the files also specify exobase, tropopause, etc.
    =#
    # do the DH altitudinal profile plot
    # set up plot
    fig, ax = subplots(figsize=(6,4))
    grid(zorder=0, color="gainsboro", which="both")
    ax2 = ax.twiny()
    subplots_adjust(wspace=0, bottom=0.15)
    ax.set_xlabel("D/H ratio (in atomic D, H)", color="teal")
    ax.set_ylabel("Altitude (km)")
    ax.set_yticks(ticks=collect(0:50:200))
    ax.tick_params(axis="x", labelcolor="teal", which="both")
    ax2.set_xlabel("HD mixing ratio", color="orange")
    ax2.tick_params(axis="x", labelcolor="orange", which="both")
    
    # Set up line colors
    colors_dumb = [cgrad(:viridis)[x] for x in range(0, stop=1, length=length(experiments))]
    colors_dumb2 = [cgrad(:plasma)[x] for x in range(0, stop=1, length=length(experiments))]
    c = Array{Float64}(undef, length(experiments), 3)
    c2 = Array{Float64}(undef, length(experiments), 3)
    # Julia generates a stupid object. convert into non-stupid object for matplotlib
    for i in range(1, length=length(colors_dumb))
        c[i, 1] = red(colors_dumb[i])
        c[i, 2] = green(colors_dumb[i])
        c[i, 3] = blue(colors_dumb[i])
        c2[i, 1] = red(colors_dumb2[i])
        c2[i, 2] = green(colors_dumb2[i])
        c2[i, 3] = blue(colors_dumb2[i])
    end

    # do actual plotting
    for i in range(1, length=length(experiments))
        ax.plot(DHproflist[i, :], alt[2:end-1]./1e5, zorder=10, color=c[i, :], 
                linewidth=3, label=optlegend*"="*experiments[i])
        ax2.plot(HDproflist[i, :], alt[2:end-1]./1e5, zorder=10, color=c2[i, :], 
                linewidth=3, label=optlegend*"="*experiments[i])
    end

    # set savepath
    plotpath = "../Results/TradeoffPlots/"*varyable*"_tradeoff_plots/"
    savepath = plotpath*varyable*optext*"_DH_and_HD_prof.png"
    #legend(fontsize=12, bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
    savefig(savepath, bbox_inches="tight")

    # save it again but with log x-axis 
    xscale("log")
    savepath = plotpath*varyable*optext*"_DH_and_HD_prof_LOG.png"
    savefig(savepath, bbox_inches="tight")
    close(fig)
end

function DH_alt_prof_plot(DHproflist, exps, v, s, optext="", optlegend="")
    #=
    DHproflist: An array with D/H profiles in each row
    exps: a list of experiment identifiers as strings
    v: water, Oflux, or temp (for putting in correct folder)
    s: a subfolder "abs/" or "mr/" for either absolute abundances or mixing ratio
    optext: an optional extension that will be appended to v, for 
            the temperature case where it goes in the temp_tradeoff_plots 
            folder but the files also specify exobase, tropopause, etc.
    optlegend: an optional string to put into the legend
    =#

    # do the DH altitudinal profile plot
    # set up plot
    fig, ax = subplots(figsize=(6,4))
    grid(zorder=0, color="gainsboro", which="both")
    subplots_adjust(wspace=0, bottom=0.15)
    ax.set_xlabel("D/H ratio (in atomic D, H)", color="teal")
    ax.set_ylabel("Altitude (km)")
    ax.set_yticks(ticks=collect(0:50:200))
    ax.tick_params(axis="x", labelcolor="teal", which="both")
    
    # generate colors
    c = get_colors(length(exps))

    # do actual plotting
    for i in range(1, length=length(exps))
        ax.plot(DHproflist[i, :], alt[2:end-1]./1e5, zorder=10, color=c[i, :], 
                linewidth=3, label=optlegend*"="*exps[i])
    end

    # set savepath
    plotpath = "../Results/TradeoffPlots/"*v*"_tradeoff_plots/"*s
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

    # where to plot the extra values 

    # Pretty ugly, but data structs are always better than if/else, right?
    paststudies = Dict("water"=>Dict("yung"=>7, "nair"=>3),
                            "Oflux"=>Dict("yung"=>1, "nair"=>5),
                            "temp_surface"=>Dict("yung"=>220, "nair"=>214), 
                            "temp_tropopause"=>Dict("yung"=>140, "nair"=>140), 
                            "temp_exobase"=>Dict("yung"=>364, "nair"=>288))

    lookupkey = pathkey*tempkey
    ystr = s == "abs/" ? L"Abundance, CO and O$_2$" : L"Mixing ratio, CO and O$_2$"
    
    # CO, O2, CO/O2 plot
    # set up plot
    fig, ax = subplots(figsize=(6,4))
    #grid(zorder=0, color="gainsboro", which="both")
    xlabel(xlab)
    xticks(rotation=45)
    ax2 = ax.twinx() # set up the other axis

    # CO and O2 mixin ratio axis
    ax.plot(xvals, ydict["CO"], marker="o", zorder=10, color="cornflowerblue") 
    ax.plot(xvals, ydict["O2"], marker="o", zorder=10, color="navy")
    ax.tick_params(axis="y", labelcolor="navy", which="both")
    ax.axvline(meanX, color="black", label="Nominal value")
    ax.set_ylabel(ystr, color="navy")
    ax.set_yscale("log")
    
    # CO/O2 axis
    ax2.plot(xvals, ydict["CO/O2"], marker="o", zorder=10, color="red")
    ax2.tick_params(axis="y", labelcolor="red")
    ax2.set_ylabel("CO/O"*L"_2"*" ratio", color="red")
    ytix = collect(0:0.2:1.2)
    if maximum(ydict["CO/O2"]) / minimum(ydict["CO/O2"]) >= 10
        ax2.set_yscale("log")
    end

    println("at past studies...")
    println(lookupkey)
    # past studies
    ax2.errorbar(paststudies[lookupkey]["nair"], 0.14, capsize=2, zorder=10,
                 yerr=reshape([0.05; 0.43], 2, 1), color="xkcd:brick red",
                 ecolor="xkcd:brick red", marker="s")  # Nair 1994 low water
    ax2.scatter(paststudies[lookupkey]["yung"], 0.53, marker="*", s=100,
                color="xkcd:tangerine", zorder=10)  # Yung 1988

    # more ticks for exobase surface temp
    if lookupkey == "temp_exobase"
        ax.set_xticks(ticks=vcat(xvals, [325, 350]))
        [l.set_visible(false) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 2 != 0]
    end

    plotpath = "../Results/TradeoffPlots/"*pathkey*"_tradeoff_plots/"*s
    savepath = plotpath*lookupkey*"_CO_and_O2.png"
    savefig(savepath, bbox_inches="tight")
    close(fig)
end

function get_colors(L)
    #=
    Generates some colors based on a color map for use in plotting a bunch of
    lines all at once.
    L: number of colors to generate.
    =#

    colors_dumb = [cgrad(:viridis)[x] for x in range(0, stop=1, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = red(colors_dumb[i])
        c[i, 2] = green(colors_dumb[i])
        c[i, 3] = blue(colors_dumb[i])
    end
    return c
end

function rel_change_plot(d, x, xlab, folderpath, FNext="", ylog="yes")
    #=
    d: a dictionary, keys are particular metrics varied, values are lists of the
       values of that metric for different variations of whatever x is (water, 
       O flux etc), but normalized to the value in the nominal case.
    x: x data, probably strings. maybe not though
    xlab: a label for the x axis
    folderpath: should be in the form parent/child, i.e. water_tradeoff_plots/abs.
    =#
    fig = figure(figsize=(6,4))
    grid(zorder=0, color="gainsboro")
    c = get_colors(length(keys(d)))
    marks = [".","o","v","^","<",">","1","2","3","4","8","s","p","P","*","h",
             "H","+","x","X","D","d","|","_"]

    i = 1
    for k in keys(d)
        plot(x, d[k], zorder=10, color=c[i, :], label=k, marker=marks[i])
        i +=1 
    end
    axhline(1, color="black")
    xlabel(xlab)
    xticks(rotation=45)
    ylabel("Rel. change vs nominal")
    if ylog=="yes"
        yscale("log")
    end
    legend(bbox_to_anchor=(1.1, 1))
    plotpath = "../Results/TradeoffPlots/"*folderpath
    savepath = "relchange"*FNext*".png"
    savefig(plotpath*savepath, bbox_inches="tight")
    close(fig)
end

# fundamental constants ========================================================
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm
DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Basic altitude stuff
alt = (0:2e5:200e5)
n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
const zmax = alt[end];

# Functions to make the plots ==================================================

# foolder location
# lead = "/data/GDrive-CU/Research/Results/TradeoffPlots/"
lead = "/home/emc/GDrive-CU/Research/Results/TradeoffPlots/"

function make_water_var_plots(abs_or_mr, path=lead)
    # Establish parameters, filenames, etc
    watervals_str = ["1e-6", "5e-6", "1e-5", "5e-5", "1e-4", "5e-4", "1e-3", "5e-3"]
    wfilelist = [path*"water_"*w*"/converged_standardwater_D_water_"*w*".h5" for w in watervals_str]
    temps = [192.0, 110.0, 199.0]
    oflux = 1.2e8
    nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
    q = abs_or_mr == "abs" ? " abundance" : " mixing ratio" # for labels
    mean_idx = findfirst(isequal("1e-4"), watervals_str) - 1
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"
    
    # Establish variables to store data on simulations
    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(watervals_str), length(alt)-2)  
    HDprofs = Array{Any}(undef, length(watervals_str), length(alt)-2)
    wdict = Dict{String, Array}("O2"=>[], "HD"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                 "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                 "DH"=>[])
    normed_dict = Dict()
    individual_plots = filter!(e->e∉["CO", "O2", "CO/O2"], [k for k in keys(wdict)])

    # loop water files and collect data 
    i = 1
    for wfile in wfilelist
        # get the current array
        ncur = get_ncurrent(wfile)

        N = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)

        # Calculate the things we care about
        # O2 Mixing ratio at surface
        append!(wdict["O2"], ncur[:O2][1]/N)
        # HD mixing ratio
        append!(wdict["HD"], ncur[:HD][1]/N)
        # H2 mixing ratio
        append!(wdict["H2"], ncur[:H2][1]/N)
        # D Mixing ratio
        append!(wdict["D"], ncur[:D][end]/N)
        # H mixing ratio
        append!(wdict["H"], ncur[:H][end]/N)
        # D/H at 150 km
        append!(wdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # CO mixing ratio
        append!(wdict["CO"], ncur[:CO][1]/N)
        # CO/O2 ratio
        append!(wdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 mixing ratio
        append!(wdict["O3"], ncur[:O3][1]/N)
        # H and D fluxes
        Hf = get_H_fluxes(wfile, oflux, temps)
        Df = get_D_fluxes(wfile, oflux, temps)
        append!(wdict["Hflux"], Hf)
        append!(wdict["Dflux"], Df)
        # fractionation factor 
        append!(wdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
        HDprofs[i, :] = ncur[:HD] ./ [n_tot(ncur, i) for i in alt[2:end-1]]
        i += 1
    end

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Laksaman"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = "Water"*q*", lower atmo"

    # Do the relative increase plot
    for k in keys(wdict)
        base_i = findfirst(isequal("1e-4"), watervals_str)
        normed_dict[k] = normalize(wdict[k], base_i)
    end

    if abs_or_mr == "mr"
        rel_change_plot(normed_dict, watervals_str, xlab, "water_tradeoff_plots/")
    end

    # Do the individual plots
    for i in individual_plots
        # set up plot
        fig, ax = subplots(figsize=(6,4))
        grid(zorder=0, color="gainsboro")
        plot(watervals_str, wdict[i], marker="o", zorder=10) 
        ax.axvline(mean_idx, color="black", label="Nominal value")

        # past studies
        if i=="f"
            scatter(findfirst(isequal("5e-3"), watervals_str)-1, 0.32, marker="*", 
                    s=100, color="green", zorder=10)  # Yung 1988
        end
        
        # set axes labels
        xlabel("Water mixing ratio, lower atmo")
        xticks(rotation=45)
        ylabdict = Dict("DH"=>"D/H ratio (/1.6e-4) @ 150 km",
                        "D"=>"D"*q*" at exobase",
                        "H"=>"H"*q*" at exobase",
                        "Dflux"=>L"$\phi_D$ (cm$^{-2}$s$^{-1}$)",
                        "Hflux"=>L"$\phi_H$ (cm$^{-2}$s$^{-1}$)",
                        "f"=>"fractionation factor (f)",
                        "O2"=>"O"*L"_2"*q,
                        "HD"=>"HD"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>"O"*L"_3"*q)
        ylabel(ylabdict[i])

        # manual ylim and yscale for things not in nologplease
        if ~(i in nologplease)#maximum(odict[i])/minimum(odict[i]) > 10#
            yscale("log")
            if i in ["Hflux", "HD", "H", "H2", "O3", "O2"]
                ylim(minimum(wdict[i])/2, maximum(wdict[i])*2)
            end
        end

        # set savepath
        plotpath = "../Results/TradeoffPlots/water_tradeoff_plots/" * subfolder
        savepath = "water_"*i*".png"
        
        savefig(plotpath*savepath, bbox_inches="tight")
        close(fig)
    end

    # CO, O2, and CO/O2 plot
    CO_O2_plot(watervals_str, wdict, xlab, "water", mean_idx, subfolder)

    # DH altitude plot
    DH_alt_prof_plot(DHprofs, watervals_str, "water", subfolder)

    println("Finished water plots")
end


function make_Oflux_var_plots(abs_or_mr, path=lead)
    # Establish important parameters, files, etc
    Ofluxvals = [8e7, 9e7, 1e8, 1.1e8, 1.2e8, 1.3e8, 1.4e8, 1.5e8, 1.6e8]
    Ofluxvals_str = ["8e7", "9e7", "1e8", "1.1e8", "1.2e8", "1.3e8", "1.4e8", 
                     "1.5e8", "1.6e8"]
    Ofilelist = [path*"Oflux_"*o*"/converged_standardwater_D_Oflux_"*o*".h5" 
                 for o in Ofluxvals_str]
    temps = [192.0, 110.0, 199.0]
    nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
    mean_idx = findfirst(isequal(1.2e8), Ofluxvals) - 1
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # Establish variables to store data on simulations
    odict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                 "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                 "DH"=>[])
    normed_dict = Dict()
    individual_plots = filter!(e->e∉["CO", "O2", "CO/O2", "DH", "HD"], 
                                  [k for k in keys(odict)])

    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(Ofluxvals_str), length(alt)-2)
    #HDprofs = Array{Any}(undef, length(Ofluxvals_str), length(alt)-2)

    # loop through O flux files 
    i = 1
    for (oflux, ofile) in zip(Ofluxvals, Ofilelist)
        # calculate the things we care about
        # get the current array
        ncur = get_ncurrent(ofile)

        N = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)

        # Calculate the things we care about
        # O2 Mixing ratio at surface
        append!(odict["O2"], ncur[:O2][1]/N)
        # HD mixing ratio
        append!(odict["HD"], ncur[:HD][1]/N)
        # H2 mixing ratio
        append!(odict["H2"], ncur[:H2][1]/N)
        # D Mixing ratio
        append!(odict["D"], ncur[:D][1]/N)
        # H mixing ratio
        append!(odict["H"], ncur[:H][1]/N)
        # D/H at 150 km
        append!(odict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # CO mixing ratio
        append!(odict["CO"], ncur[:CO][1]/N)
        # CO/O2 ratio
        append!(odict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 mixing ratio
        append!(odict["O3"], ncur[:O3][1]/N)
        # H and D fluxes
        Hf = get_H_fluxes(ofile, oflux, temps)
        Df = get_D_fluxes(ofile, oflux, temps)
        append!(odict["Hflux"], Hf)
        append!(odict["Dflux"], Df)
        # fractionation factor 
        append!(odict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profiles
        #HDprofs[i, :] = ncur[:HD] ./ [n_tot(ncur, i) for i in alt[2:end-1]]
        i += 1
    end

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Laksaman"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    xlab = L"$\phi_O$ (cm$^{-2}$s$^{-1}$)"

    # Do the relative increase plot
    for k in keys(odict)
        base_i = findfirst(isequal(1.2e8), Ofluxvals)
        normed_dict[k] = normalize(odict[k], base_i)
    end

    if abs_or_mr == "mr"
        rel_change_plot(normed_dict, Ofluxvals_str, xlab, "Oflux_tradeoff_plots/", "no")
    end

    # Make individual plots showing change 
    for i in individual_plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        grid(zorder=0, color="gainsboro")
        plot(Ofluxvals_str, odict[i], marker="o", zorder=10) 
        ax.axvline(mean_idx, color="black", label="Nominal value")

        # past studies
        if i=="f"
            scatter(findfirst(isequal(8e7), Ofluxvals)-1, 0.32, marker="*", s=100,
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

        # manual ylim and yscale for things not in nologplease
        if ~(i in nologplease)#maximum(odict[i])/minimum(odict[i]) > 10#
            yscale("log")
            if i in ["Hflux", "HD", "O3"]
                ylim(minimum(odict[i])/2, maximum(odict[i])*2)
            end
        end

        # set savepath
        plotpath = "../Results/TradeoffPlots/Oflux_tradeoff_plots/"*subfolder
        savepath = plotpath*"O_flux_"*i*".png"
        
        savefig(savepath, bbox_inches="tight")
        close(fig)
    end

    # CO, O2 and CO/O2 Plot
    CO_O2_plot(Ofluxvals_str, odict, xlab, "Oflux", 
               mean_idx, subfolder)

    DH_alt_prof_plot(DHprofs, Ofluxvals_str, "Oflux", subfolder)

    println("Finished O flux")
end
    

function make_T_var_plots(abs_or_mr, path=lead)

    # Set up parameters, filenames, etc
    tvals = Dict("surface"=>[180.0, 190.0, 200.0, 210.0, 
                             220.0, 230.0, 240.0, 250.0, 260.0, 270.0],
                 "tropopause"=>[70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 
                                140.0, 150.0, 160.0],
                 "exobase"=>[125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 
                             275.0, 300.0])
    tvals_str = Dict()
    tempfilelist = Dict()
    for k in keys(tvals)
        tvals_str[k] = [string(trunc(Int, x)) for x in tvals[k]]
        if k == "surface"
            tempfilelist[k] = [path*"temp_"*t*"_110_199"*"/converged_standardwater_D_temp_"*t*"_110_199.h5" for t in tvals_str[k]]
        elseif k == "tropopause"
            tempfilelist[k] = [path*"temp_192_"*t*"_199"*"/converged_standardwater_D_temp_192_"*t*"_199.h5" for t in tvals_str[k]]
        elseif k == "exobase"
            tempfilelist[k] = [path*"temp_192_110_"*t*"/converged_standardwater_D_temp_192_110_"*t*".h5" for t in tvals_str[k]]
        end
    end
    meanT = Dict("surface"=>190, "tropopause"=>110, "exobase"=>200)  # nominal 
    oflux = 1.2e8
    nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # loop through which temp is varied and construct a list of datapoints
    for experiment in keys(tvals) # loop across the dictionary
        println("Now doing temperature ", experiment)
        tdict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                     "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                     "DH"=>[])
        normed_dict = Dict()
        DHprofs = Array{Any}(undef, length(tvals_str[experiment]), length(alt)-2)
        #HDprofs = Array{Any}(undef, length(tvals_str[experiment]), length(alt)-2)
        individual_plots = filter!(e->e∉["CO", "O2", "CO/O2", "DH"], 
                                        [k for k in keys(tdict)])

        # now loop through the values for each varied temp
        i = 1
        for (tv, tfile) in zip(tvals[experiment], tempfilelist[experiment])
            # set the temperature profile
            if experiment == "surface" 
                temps = [tv, 110.0, 199.0]
            elseif experiment == "tropopause"
                temps = [192.0, tv, 199.0]
            elseif experiment == "exobase"
                temps = [192.0, 110.0, tv]
            end

            # get the current array
            ncur = get_ncurrent(tfile)

            N = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0)

            # Calculate the things we care about
            # O2 Mixing ratio at surface
            append!(tdict["O2"], ncur[:O2][1]/N)
            # HD mixing ratio
            append!(tdict["HD"], ncur[:HD][1]/N)
            # H2 mixing ratio
            append!(tdict["H2"], ncur[:H2][1]/N)
            # D Mixing ratio
            append!(tdict["D"], ncur[:D][1]/N)
            # H mixing ratio
            append!(tdict["H"], ncur[:H][1]/N)
            # D/H at 150 km
            append!(tdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
            # CO mixing ratio
            append!(tdict["CO"], ncur[:CO][1]/N)
            # CO/O2 ratio
            append!(tdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
            # O3 mixing ratio
            append!(tdict["O3"], ncur[:O3][1]/N)
            # H and D fluxes
            Hf = get_H_fluxes(tfile, oflux, temps)
            Df = get_D_fluxes(tfile, oflux, temps)
            append!(tdict["Hflux"], Hf)
            append!(tdict["Dflux"], Df)
            # fractionation factor 
            append!(tdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
            # D/H profile
            DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
            #HDprofs[i, :] = ncur[:HD] ./ [n_tot(ncur, i) for i in alt[2:end-1]]
            i += 1
        end

        # make plots pretty
        rcParams = PyDict(matplotlib."rcParams")
        rcParams["font.sans-serif"] = ["Laksaman"]
        rcParams["font.monospace"] = ["FreeMono"]
        rcParams["font.size"] = 22
        rcParams["axes.labelsize"]= 24
        rcParams["xtick.labelsize"] = 22
        rcParams["ytick.labelsize"] = 22
        
        xlab = experiment*" temperature"

        # Do the relative increase plot
        for k in keys(tdict)
            base_i = findfirst(isequal(meanT[experiment]), tvals[experiment])
            normed_dict[k] = normalize(tdict[k], base_i)
        end

        if abs_or_mr == "mr"
            rel_change_plot(normed_dict, tvals_str[experiment], xlab, 
                            "temp_tradeoff_plots/", "_temp_"*experiment)# "no")
        end

        # loop through the parameters of interest and plot them
        for i in individual_plots

            # basic plot stuff
            fig, ax = subplots(figsize=(6,4))
            grid(zorder=0, color="gainsboro")
            plot(tvals[experiment], tdict[i], marker="o", zorder=10) 
            ax.axvline(meanT[experiment], color="black", label="Nominal value")
            xlabel(xlab)
            xticks(ticks=tvals[experiment], labels=tvals_str[experiment], 
                   rotation=45)
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
            if experiment=="surface"
                if i == "f"
                    scatter(220, 0.32, marker="d", color="xkcd:tangerine", 
                            zorder=10)  # Yung 1988
                end
            elseif experiment=="tropopause"
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
            elseif experiment=="exobase"
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
                    #xlim(100,400)
                    xticks(ticks=vcat(tvals[experiment], [325, 350]), 
                           labels=vcat(tvals_str[experiment], 
                           ["325", "350"]), rotation=45)
                end
            end 

            # manual ylim and yscale for things not in nologplease
            if ~(i in nologplease)#maximum(tdict[i])/minimum(tdict[i]) > 10#
                yscale("log")
                if i in ["Hflux", "HD"]
                    ylim(minimum(tdict[i])/2, maximum(tdict[i])*2)
                end
            end

            # set savepath
            plotpath = "../Results/TradeoffPlots/temp_tradeoff_plots/"*subfolder
            savepath = plotpath*join(["temp", experiment], "_")*"_"*i*".png"
            savefig(savepath, bbox_inches="tight")
            close(fig)
        end

        # CO/O2 plot
        CO_O2_plot(tvals[experiment], tdict, "Temperature at "*experiment, 
                   "temp", meanT[experiment], subfolder, "_"*experiment)

        # D/H altitude plot
        DH_alt_prof_plot(DHprofs, tvals_str[experiment], "temp", subfolder, "_"*experiment, 
                         latexstring("T_{$(experiment[1:3])}"))

        println("Finished temperature ", experiment)
    end
    println("Finished temps")
end

# make_water_var_plots("mr")
# make_water_var_plots("abs")
# println()
# make_Oflux_var_plots("abs")
# make_Oflux_var_plots("mr")
# println()
make_T_var_plots("mr")