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
using Analysis
using DataFrames

include("PARAMETERS.jl")

# fundamental constants ========================================================
DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

global tvals = Dict("surface"=>[150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 
                                230.0, 240.0, 250.0, 260.0, 270.0],
                    "tropopause"=>[100.0, 110.0, 120.0, 130.0, #70.0, 80.0, 90.0, 
                                   140.0, 150.0, 160.0],
                    "exobase"=>[150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 
                                300.0, 325.0, 350.0])
global Oflux_vals = ["8e7", "9e7", "1.0e8", "1.1e8", "1.2e8", "1.3e8", 
                            "1.4e8", "1.5e8", "1.6e8"]

global watervals = [1, 10, 25, 50, 100]#, 150]
global watervals_str = ["1.33e-5", "1.41e-4", "3.72e-4", "8.2e-4", "2.89e-3"]#, "1.358e-2"]

# nominal value plot location or index, 1-indexed for Julia
global nom_i_julia = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs, 
                   "O flux"=>findfirst(isequal("1.2e8"), Oflux_vals), 
                   "water"=>10)
# Passing the values to PyPlot requires 0-indexing so here it is:
global nom_i_py = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs, 
                "O flux"=>findfirst(isequal("1.2e8"), Oflux_vals)-1, 
                "water"=>10)

# set the data for comparison. Order:
# CO MR (Trainer+2019), O2 MR (Trainer+2019), H2 abundance (Kras&Feldman 2001),
# O3 μm-atm (Clancy 2016)
global data = [5.8e-4, 1.61e-3, 15, 1.18]  # CO MR, O2 MR, H2 ppm, O3 μm-atm
global s = [0.8e-4, 0.09e-3, 5, 0.7]     # sigmas (uncertainty) on each

medgray = "#444444"
sz = 10

# Functions ====================================================================

function normalize(arr, base_i)
    normed_arr = arr ./ arr[base_i]
    return normed_arr
end

function normalize_val(arr, normval)
    normed_arr = arr ./ normval
    return normed_arr
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
    plot_bg(ax)
    subplots_adjust(wspace=0, bottom=0.15)
    ax.set_xlabel("D/H ratio (in atomic D, H)")
    ax.set_ylabel("Altitude (km)")
    ax.set_yticks(ticks=collect(0:50:200))
    
    # generate colors
    c = get_grad_colors(length(exps), "viridis")

    # do actual plotting
    if typeof(exps[1]) != String
        exps = [string(x) for x in exps]
    end
    for i in range(1, length=length(exps))
        ax.plot(DHproflist[i, :], alt[2:end-1]./1e5, zorder=10, color=c[i, :], 
                linewidth=3, label=optlegend*"="*exps[i])
    end

    # set savepath
    plotpath = mainpath*v*"_plots/"*s
    savepath = plotpath*v*optext*"_DH_prof.png"
    legend(fontsize=12, bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
    savefig(savepath, bbox_inches="tight")

    # save it again but with log x-axis 
    xscale("log")
    # xlim(3.5e-4,5e-4)
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
    paststudies = Dict("water"=>Dict("yung"=>15, "nair"=>[3,8.8]), # these are in pr μm
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

    plotpath = mainpath*pathkey*"_plots/"*s
    savepath = plotpath*lookupkey*"_CO_and_O2.png"
    savefig(savepath, bbox_inches="tight")
    close(fig)
end

# For the main case in temperatures since they have weird numbers --------------

function makenomdict(abs_or_mr)
    #=
    TODO: check this isn't broken
    =#

    # get the current array
    tfile = mainpath*"temp_$(meanTsint)_$(meanTtint)_$(meanTeint)/converged_temp_$(meanTsint)_$(meanTtint)_$(meanTeint).h5"
    ncur =  get_ncurrent(tfile)
    N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0, n_alt_index)
    Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, zmax, n_alt_index)
    LA = collect(0e5:2e5:78e5)

    nomTdict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H2MR"=>[], "H"=>[], "D"=>[], "CO"=>[], 
                     "Hflux"=>[], "Dflux"=>[], "f"=>[], "CO/O2"=>[], "O3"=>[], 
                     "DH"=>[])
    DHprofs = Array{Any}(undef, 1, length(alt)-2)
    # Calculate the things we care about
    # O2 Mixing ratio at surface
    append!(nomTdict["O2"], ncur[:O2][1]/N0)
    # HD mixing ratio
    append!(nomTdict["HD"], ncur[:HD][end]/Ntop)
    # H2 mixing ratio
    append!(nomTdict["H2MR"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h, n_alt_index) for h in LA]))
    append!(nomTdict["H2"], ncur[:H2][end]/Ntop)
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
    Hf = get_flux(:H, "thermal", tfile, 1.2e8, meantemps)
    Df = get_flux(:D, "thermal", tfile, 1.2e8, meantemps)
    append!(nomTdict["Hflux"], Hf)
    append!(nomTdict["Dflux"], Df)
    # fractionation factor 
    append!(nomTdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
    # D/H profile
    DHprofs[1, :] = ncur[:D] ./ ncur[:H]  # altitude profile

    return nomTdict
end

# Analyzation functions (main routines) ----------------------------------------

function analyze_water(abs_or_mr, allDbearers, make_plots=false, path=mainpath)
    # Establish parameters, filenames, etc
    wfilelist = [path*"water_"*w*"/converged_water_"*w*".h5" for w in watervals_str]
    temps = [meanTs, meanTt, meanTe]
    oflux = 1.2e8
    q = abs_or_mr == "abs" ? " abundance" : " mixing ratio" # for labels
    mean_idx = findfirst(isequal(10), watervals) - 1 
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"
    
    # Establish variables to store data on simulations
    # Easier to deal with D/H profiles separately due to different array size
    DHprofs = Array{Any}(undef, length(watervals), length(alt)-2)  
    wdict = Dict{String, Array}("O2"=>[], "HD"=>[], "H2"=>[], "H2MR"=>[], "H"=>[], "D"=>[], "CO"=>[], 
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

        N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0, n_alt_index)
        Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, 250e5, n_alt_index)
        LA = collect(0e5:2e5:78e5)

        # Calculate the things we care about
        # O2 Mixing ratio at surface
        append!(wdict["O2"], ncur[:O2][1]/N0)
        # HD at exobase
        append!(wdict["HD"], ncur[:HD][end]/Ntop)
        # H2 mixing ratio in lower atmo and at exobase
        append!(wdict["H2MR"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h, n_alt_index) for h in LA])) # ppm in lower atmo
        append!(wdict["H2"], ncur[:H2][end]/Ntop)
        # D at exobase
        append!(wdict["D"], ncur[:D][end]/Ntop)
        # H at exobase
        append!(wdict["H"], ncur[:H][end]/Ntop)
        # D/H at 150 km
        append!(wdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # CO mixing ratio at surface 
        append!(wdict["CO"], ncur[:CO][1]/N0)
        # CO/O2 ratio at surface
        append!(wdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 in #/cm^2, used to convert to μm-atm later
        append!(wdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
        # H and D fluxes
        Hf = get_flux(:H, "thermal", wfile, oflux, temps)
        Df = get_flux(:D, "thermal", wfile, oflux, temps)
        append!(wdict["Hflux"], Hf)
        append!(wdict["Dflux"], Df)
        # fractionation factor 
        append!(wdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
        # Other D bearing species at exobase
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

    return wdict
end

function analyze_Oflux(abs_or_mr, allDbearers, make_plots=false, path=mainpath)
    # Establish important parameters, files, etc
    Ofluxvals = [8e7, 9e7, 1e8, 1.1e8, 1.2e8, 1.3e8, 1.4e8, 1.5e8, 1.6e8]
    Ofluxvals_str = ["8e7", "9e7", "1.0e8", "1.1e8", "1.2e8", "1.3e8", "1.4e8", "1.5e8", "1.6e8"]
    Ofilelist = [path*"Oflux_"*o*"/converged_Oflux_"*o*".h5" 
                 for o in Ofluxvals_str]
    temps = [meanTs, meanTt, meanTe]
    mean_idx = findfirst(isequal(1.2e8), Ofluxvals) - 1
    q = abs_or_mr == "abs" ? " abundance " : " mixing ratio " # set label
    subfolder = abs_or_mr == "abs" ? "abs/" : "mr/"

    # Establish variables to store data on simulations
    odict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H2MR"=>[], "H"=>[], "D"=>[], "CO"=>[], 
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
        # get the current array
        ncur = get_ncurrent(ofile)

        N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0, n_alt_index)
        Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, zmax, n_alt_index)
        LA = collect(0e5:2e5:78e5)

        # Calculate the things we care about
        # O2 Mixing ratio at surface
        append!(odict["O2"], ncur[:O2][1]/N0)
        # HD at exobase
        append!(odict["HD"], ncur[:HD][end]/Ntop)
        # H2 mixing ratio in lower atmo and at exobase
        append!(odict["H2MR"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h, n_alt_index) for h in LA])) # ppm in lower atmo
        append!(odict["H2"], ncur[:H2][end]/Ntop)
        # D at exobase
        append!(odict["D"], ncur[:D][end]/Ntop)
        # H at exobase
        append!(odict["H"], ncur[:H][end]/Ntop)
        # D/H at 150 km
        append!(odict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
        # CO mixing ratio at surface 
        append!(odict["CO"], ncur[:CO][1]/N0)
        # CO/O2 ratio at surface
        append!(odict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
        # O3 in #/cm^2, used to convert to μm-atm later
        append!(odict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
        # H and D fluxes
        Hf = get_flux(:H, "thermal", ofile, oflux, temps)
        Df = get_flux(:D, "thermal", ofile, oflux, temps)
        append!(odict["Hflux"], Hf)
        append!(odict["Dflux"], Df)
        # fractionation factor 
        append!(odict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
        # D/H profile
        DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
        # Other D bearing species abundances at exobase
        if allDbearers
            append!(odict["OD"], ncur[:OD][end]/Ntop)
            append!(odict["HDO2"], ncur[:HDO2][end]/Ntop)
            append!(odict["DO2"], ncur[:DO2][end]/Ntop)
        end
        i += 1
    end

    if make_plots == true
        make_Oflux_plots(Ofluxvals, Ofluxvals_str, odict, DHprofs, q, mean_idx, 
                         subfolder)
    end

    return odict
end
    
function analyze_T(abs_or_mr, allDbearers, make_plots=false, path=mainpath)
    #=
    abs_or_mr: whether using the absolute abundance file or mixing ratio file
    allDbearers: whether to plot extra D-bearing species
    make_plots: whether to generate the plots
    path: where to save plots
    =#

    # Set up parameters, filenames, etc
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
        tdict = Dict("O2"=>[], "HD"=>[], "H2"=>[], "H2MR"=>[], "H"=>[], "D"=>[], "CO"=>[], 
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
            N0 = abs_or_mr == "abs" ? 1 : n_tot(ncur, 0, n_alt_index)
            Ntop = abs_or_mr == "abs" ? 1 : n_tot(ncur, zmax, n_alt_index)
            LA = collect(0e5:2e5:78e5)

            # Calculate the things we care about
            # O2 Mixing ratio at surface
            append!(tdict["O2"], ncur[:O2][1]/N0)
            # HD at exobase
            append!(tdict["HD"], ncur[:HD][end]/Ntop)
            # H2 mixing ratio in lower atmo and at exobase
            append!(tdict["H2MR"], sum(ncur[:H2][1:length(LA)])/sum([n_tot(ncur, h, n_alt_index) for h in LA])) # ppm in lower atmo
            append!(tdict["H2"], ncur[:H2][end]/Ntop)
            # D at exobase
            append!(tdict["D"], ncur[:D][end]/Ntop)
            # H at exobase
            append!(tdict["H"], ncur[:H][end]/Ntop)
            # D/H at 150 km
            append!(tdict["DH"], (ncur[:D][75]/ncur[:H][75])/(1.6e-4))
            # CO mixing ratio at surface 
            append!(tdict["CO"], ncur[:CO][1]/N0)
            # CO/O2 ratio at surface
            append!(tdict["CO/O2"], ncur[:CO][1]/ncur[:O2][1])
            # O3 in #/cm^2, used to convert to μm-atm later
            append!(tdict["O3"], sum(ncur[:O3])*2e5) # gets O3 in #/cm^2
            # H and D fluxes
            Hf = get_flux(:H, "thermal", tfile, oflux, temps)
            Df = get_flux(:D, "thermal", tfile, oflux, temps)
            append!(tdict["Hflux"], Hf)
            append!(tdict["Dflux"], Df)
            # fractionation factor 
            append!(tdict["f"], 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1]))
            # D/H profile
            DHprofs[i, :] = ncur[:D] ./ ncur[:H]  # altitude profile
            # Other D bearing species at exobase
            if allDbearers
                append!(tdict["OD"], ncur[:OD][end]/Ntop)
                append!(tdict["HDO2"], ncur[:HDO2][end]/Ntop)
                append!(tdict["DO2"], ncur[:DO2][end]/Ntop)
            end
            i += 1
        end

        if make_plots == true
            make_T_plots(tvals, tvals_str, tdict, DHprofs, experiment, q, meanT, 
                         subfolder)
        end

        all_tdicts[experiment] = tdict
    end
    return all_tdicts
end

# Main specific plots of various species with each parameter -------------------
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

    xlab = "Water"*q*L" (pr $\mu$m)"

    # Do the individual plots
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HD", "OD", "DO2", "HDO2"], [k for k in keys(d)])
    for i in plots
        # set up plot
        fig, ax = subplots(figsize=(6,4))
        plot_bg(ax)
        xscale("log")
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
                        "H2MR"=>"H"*L"_2"*q,
                        "H2"=>"H"*L"_2"*q,
                        "O3"=>L"O$_3$ (#/cm$^{-2}$)")
        ylabel(ylabdict[i])

        # Certain measureables have a range not suited to log scale
        nologplease = ["Dflux", "CO/O2", "DH"]  # don't logscale everything
        if ~(i in nologplease)
            yscale("log")
            if i in ["Hflux", "HD", "H", "H2", "O3", "O2"]
                ylim(minimum(d[i])/2, maximum(d[i])*2)
            end
        end

        # set savepath
        plotpath = mainpath*"water_plots/" * s
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
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HD", "OD", "DO2", "HDO2"], [k for k in keys(d)])
    for i in plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        plot_bg(ax)
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
                        "H2MR"=>"H"*L"_2"*q,
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
        plotpath = mainpath*"Oflux_plots/"*s
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
    plots = filter!(e->e∉["CO", "O2", "CO/O2", "HD", "OD", "DO2", "HDO2"], [k for k in keys(d)])
    for i in plots
        # basic plot stuff
        fig, ax = subplots(figsize=(6,4))
        plot_bg(ax)
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
                        "H2MR"=>"H"*L"_2"*q,
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
        plotpath = mainpath*"temp_plots/"*s
        savepath = plotpath*join(["temp", exp], "_")*"_"*i*".png"
        savefig(savepath, bbox_inches="tight")
        close(fig)
    end

    # make special plots
    CO_O2_plot(T[exp], d, exp*" temperature", "temp", nomT[exp], s, "_"*exp)
    DH_alt_prof_plot(DHdata, T_str[exp], "temp", s, "_"*exp, 
                     latexstring("T_{$(exp[1:3])}"))
end

# Make the key 'trade off' plots -----------------------------------------------
function make_Oflux_main_plots(output_MR, output_abs)
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

    # first plot - compare with data -------------------------------------------
    fig, ax1 = subplots(figsize=(7,6))
    plot_bg(ax1)

    c1 = ["#10007A", "#2F7965", "#e16262", "#D59A07"]  # just colors
    
    # calculate the relativeness of each point wrt data
    COdiff = (output_MR["CO"] .- data[1])/s[1]
    O2diff = (output_MR["O2"] .- data[2])/s[2]
    # H2diff = (ABSdictvar["H2"] .- data[3])/s[3]  # absolute abundance
    H2diff = (output_MR["H2MR"]./1e-6 .- data[3])/s[3] # ppm
    O3diff = (areadensity_to_micron_atm(output_abs["O3"]) .- data[4])/s[4] 


    # plot the actual stuff   
    ax1.plot(Oflux_vals, COdiff, marker="o", ms=sz, color=c1[1], zorder=10)
    ax1.plot(Oflux_vals, O2diff, marker="x", ms=sz, color=c1[2], zorder=10)
    ax1.plot(Oflux_vals, H2diff, marker="*", ms=sz, color=c1[3], zorder=10)
    ax1.plot(Oflux_vals, O3diff, marker="^", ms=sz, color=c1[4], zorder=10)
    # plot the mean values
    ax1.axvline(nom_i_py["Oflux"], color=medgray, zorder=5)
    ax1.axhline(0, color="black")

    # Text on plots 
    ax1.text(8.7, 1, "Global\nmean", color=medgray, ha="right")

    ax1.set_xlabel(L"\Phi_O"*L" (cm$^{-2}s^{-1}$)")
    ax1.set_ylabel(L"($X_{model}$-$X_{obs}$)/$\sigma$")
    savefig(mainpath*"output_vs_data_oflux.png", bbox_inches="tight")


    # second plot: H2, HD, H, D, Hflux, Dflux ----------------------------------
    fig, ax2 = subplots(figsize=(7,6))

    c2 = ["#6270d9","#95b034","#d14f58","#5ca85c","#ce6d30", "#0d186d"]

    # to do the division/normalization we need to re-find the index of the value
    # against which to normalize. To do this, look in the xvals by experiment, 
    # then index that according to whether we cut it or not. There doesn't seem 
    # to be a cleaner way to do the indexing of cut.
    normidx = findfirst(isequal("1.2e8"), Oflux_vals[1:1:length(Oflux_vals)])
    
    HDdiff = normalize(output_MR["HD"], normidx)
    H2diff = normalize(output_MR["H2"], normidx)
    Hdiff = normalize(output_MR["H"], normidx)
    Ddiff = normalize(output_MR["D"], normidx)
    Hfdiff = normalize(output_MR["Hflux"], normidx)
    Dfdiff = normalize(output_MR["Dflux"], normidx)

    # plot the actual stuff
    ax2.plot(Oflux_vals, Hdiff, marker="x", ms=sz,  color=c2[2], zorder=10, label="H")
    ax2.plot(Oflux_vals, Ddiff, marker="*", ms=sz,  color=c2[5], zorder=10, label="D")
    ax2.plot(Oflux_vals, H2diff, marker="o", ms=sz, color=c2[1], zorder=10, label="H2")
    ax2.plot(Oflux_vals, HDdiff, marker="o", ms=sz, color=c2[6], zorder=10, label="HD")
    ax2.plot(Oflux_vals, Hfdiff, marker="^", ms=sz, color=c2[4], zorder=10, label=L"\phi_H")
    ax2.plot(Oflux_vals, Dfdiff, marker="D", ms=sz, color=c2[3], zorder=10, label=L"\phi_D")
    # nominal value
    ax2.axvline(nom_i_py["Oflux"], color=medgray, zorder=5)
    ax2.axhline(1, color=medgray, zorder=5)

    # Text on plots 
    ax2.text(8.7, 1, "Global\nmean", color=medgray, ha="right")
    
    # Other species that can be plotted, or species at other locations
    # HDtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ODtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # HDO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # DO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ax2.plot(xvals[ex], HDtopdiff, marker="o", color="#7D3403", zorder=10)
    # ax2.plot(xvals[ex], ODtopdiff, marker="o", color="red", zorder=10)
    # ax2.plot(xvals[ex], HDO2topdiff, marker="o", color="blue", zorder=10)
    # ax2.plot(xvals[ex], DO2topdiff, marker="o", color="purple", zorder=10)

    ax2.set_xlabel(L"\Phi_O"*L" (cm$^{-2}s^{-1}$)")
    ax2.set_ylabel(L"X/X$_{nominal}$")
    #ax2.legend(bbox_to_anchor=(1.05, 1))
    savefig(mainpath*"compare_nominal_Oflux.png", bbox_inches="tight")

    # # third plot: f ----------------------------------------------------------
    # only for water and O flux - the others are done in a panel
    fig, ax = subplots(figsize=(7,6))
    plot_bg(ax)

    ax.tick_params(rotation=45, axis="x")
    ax.text(8.7, 0.00112, "Global\nmean", color=medgray, ha="right")


    ax.set_ylabel(L"$f$", color="black")
    ax.set_xlabel(L"\Phi_O"*L" (cm$^{-2}s^{-1}$)")
    linecol = "cornflowerblue"

    ax.plot(Oflux_vals, output_MR["f"], marker="o", color=linecol, zorder=10)
    ax.axvline(nom_i_py["Oflux"], color=medgray, zorder=5)

        
    savefig(mainpath*"f_tradeoff_Oflux.png", bbox_inches="tight")
end

# Water plots ------------------------------------------
function make_water_Hspecies_plot(output_dict, abs_or_mr)
    ############################################################
    #=
    Plot showing H, D, HD, H2, ΦH and ΦD as relative to the global mean profile
    =#

    # to do the division/normalization we need to re-find the index of the value
    # against which to normalize. To do this, look in the xvals by experiment, 
    # then index that according to whether we cut it or not. There doesn't seem 
    # to be a cleaner way to do the indexing of cut.
    normidx = findfirst(isequal(10), watervals[1:1:length(watervals)])
    
    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    # Make the plot 
    fig, ax = subplots(figsize=(7,6))
    plot_bg(ax)
    
    HDdiff = normalize(output_dict["HD"], normidx)
    H2diff = normalize(output_dict["H2"], normidx)
    Hdiff = normalize(output_dict["H"], normidx)
    Ddiff = normalize(output_dict["D"], normidx)
    Hfdiff = normalize(output_dict["Hflux"], normidx)
    Dfdiff = normalize(output_dict["Dflux"], normidx)

    # plot the actual stuff

    # colors for each line
    c = ["#95b034", "#ce6d30", "#6270d9", "#0d186d", "#5ca85c", "#d14f58"]

    ax.plot(watervals, Hdiff, marker="x", ms=sz,  color=c[1], zorder=10, label="H")
    ax.plot(watervals, Ddiff, marker="*", ms=sz+10,  color=c[2], zorder=10, 
            label="D", alpha=0.7)

    ax.plot(watervals, H2diff, marker="o", ms=sz, color=c[3], zorder=10, label="H2")
    ax.plot(watervals, HDdiff, marker="D", ms=sz, color=c[4], zorder=10, 
            label="HD", alpha=0.7)
    
    ax.plot(watervals, Hfdiff, marker="^", ms=sz, color=c[5], zorder=10, label=L"\phi_H")
    ax.plot(watervals, Dfdiff, marker="v", ms=sz, color=c[6], zorder=10, 
            label=L"\phi_D", alpha=0.7)
    # nominal value
    ax.axvline(nom_i_py["water"], color=medgray, zorder=5)
    ax.axhline(1, color=medgray, zorder=5)

    # Text on plots 
    ax.set_xscale("log")
    ax.text(11, 0.9, "Global\nmean", color=medgray, ha="left", va="top")
    ax.text(1, 0.99, "H", color=c[1], ha="left", va="top")
    ax.text(1, 0.87, "D", color=c[2], ha="left", va="top")
    ax.text(1, 1.02, L"H_2", color=c[3], ha="left", va="top")
    ax.text(2, 0.85, "HD", color=c[4], ha="left", va="top")
    ax.text(1.5, 0.995, L"$\phi_H$", color=c[5], ha="left", va="top")
    ax.text(1.5, 0.9, L"$\phi_D$", color=c[6], ha="left", va="top")
  
    # Other species that can be plotted, or species at other locations
    # HDtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ODtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # HDO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # DO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ax2.plot(xvals[ex], HDtopdiff, marker="o", color="#7D3403", zorder=10)
    # ax2.plot(xvals[ex], ODtopdiff, marker="o", color="red", zorder=10)
    # ax2.plot(xvals[ex], HDO2topdiff, marker="o", color="blue", zorder=10)
    # ax2.plot(xvals[ex], DO2topdiff, marker="o", color="purple", zorder=10)

    ax.set_xlabel(L"total atmospheric water (pr $\mu$m)")
    ax.set_ylabel(L"X/X(\overline{T}_{global})")
    savefig(mainpath*"compare_nominal_water_"*abs_or_mr*".png", bbox_inches="tight")
end

function make_water_output_vs_data(output_MR, output_abs)
    #=TODO: Fill me in=#
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22


    # first plot - compare with data -------------------------------------------
    fig, ax = subplots(figsize=(7,6))
    plot_bg(ax)

    c1 = ["#10007A", "#2F7965", "#e16262", "#D59A07"]  # just colors
    medgray = "#444444"
    
    # calculate the relativeness of each point wrt data
    COdiff = (output_MR["CO"] .- data[1])/s[1]
    O2diff = (output_MR["O2"] .- data[2])/s[2]
    # H2diff = (ABSdictvar["H2"] .- data[3])/s[3]  # absolute abundance
    H2diff = (output_MR["H2MR"]./1e-6 .- data[3])/s[3] # ppm
    O3diff = (areadensity_to_micron_atm(output_abs["O3"]) .- data[4])/s[4] 

    # plot the actual stuff   
    ax.plot(watervals, COdiff, marker="o", ms=sz, color=c1[1], zorder=10)
    ax.plot(watervals, O2diff, marker="x", ms=sz, color=c1[2], zorder=10)
    ax.plot(watervals, H2diff, marker="*", ms=sz, color=c1[3], zorder=10)
    ax.plot(watervals, O3diff, marker="^", ms=sz, color=c1[4], zorder=10)
    # plot the mean values
    ax.axvline(nom_i_py["water"], color=medgray, zorder=5)
    ax.axhline(0, color="black")

    # Text on plots 
    ax.set_xscale("log")
    ax.text(10.2, 2.5, "Global\nmean", color=medgray, ha="left", va="top")
    ax.text(1, -2.4, "CO", color=c1[1], ha="left", va="top")
    ax.text(20, -2.5, L"O$_2$", color=c1[2], ha="left", va="top")
    ax.text(1, -0.5, L"H$_2$", color=c1[3], ha="left", va="top")
    ax.text(1, 2.5, L"O$_3$", color=c1[4], ha="left", va="top")
    ax.set_xlabel(L"total atmospheric water (pr $\mu$m)")
    ax.set_ylabel(L"($X_{model}$-$X_{obs}$)/$\sigma$")
    savefig(mainpath*"output_vs_data_water.png", bbox_inches="tight")
end

function make_water_f_plot(output_MR)
    #=f as a function of water vapor=#
    fig, ax = subplots(figsize=(7,6))
    plot_bg(ax)
    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22


    ax.plot(watervals, output_MR["f"], marker="o", color="cornflowerblue", zorder=10)
    ax.axvline(nom_i_py["water"], color=medgray, zorder=5)

    ax.text(9, 0.0019, "Global\nmean", color=medgray, ha="right", va="bottom")
    ax.set_xscale("log")

    ax.set_ylabel(L"$f$", color="black")
    ax.set_xlabel(L"total atmospheric water (pr $\mu$m)")
    
    savefig(mainpath*"f_vs_water.png", bbox_inches="tight")
end

# Temp plots -------------------------------------------
function make_T_Hspecies_plot(output_dict, abs_or_mr)
    #=
    Plot showing H, D, HD, H2, ΦH and ΦD as relative to the global mean profile
    =#

    # where to place text, basically.
    exps = ["surface", "tropopause", "exobase"]
    gmean_txt_loc = Dict("surface"=>[219, 1.4], "tropopause"=>[132, 1.9], 
                         "exobase"=>[208, 100])
    # ylims = Dict("surface"=>[-10, 190], "tropopause"=>[-4, 4.2], "exobase"=>[-4, 4.2])
    xt_args = Dict("surface"=>150:20:280, "tropopause"=>100:10:160, 
                   "exobase"=>150:50:350)
    if abs_or_mr == "mr"
        linelbls = DataFrame(Exp=["surface", "tropopause", "exobase"], 
                              H=[[150, 1.8],     [100, 2.4],  [150, 10]], 
                              D=[[260, 0.5],     [100, 1.5],  [340, 0.1]],
                              H2=[[190, 2],      [105, 3.2], [150, 2.2]],
                              HD=[[150, 1.5],    [100, 2.9], [340, 0.8]],
                              fluxH=[[150, 1],  [100, 0.9], [340, 2.2]],
                              fluxD=[[150, 0.7], [100, 0.4],  [155, 0.4]])
    else
        linelbls = DataFrame(Exp=["surface", "tropopause", "exobase"], 
                              H=[[150, 1.15],   [150, 1.15],  [150, 20]], 
                              D=[[175, 0.75],   [108, 0.8],  [340, 0.3]],
                              H2=[[195, 1.4],   [100, 1.45], [150, 5]],
                              HD=[[165, 1.25],  [100, 1.15], [340, 2]],
                              fluxH=[[150, 1],  [155, 1.2], [150, 0.8]],
                              fluxD=[[160, 0.7], [100, 0.7],  [145, 0.1]])
    end

     # colors for each 
    c = ["#95b034", "#ce6d30", "#6270d9", "#0d186d", "#5ca85c", "#d14f58"]
    # to do the division/normalization we need to re-find the index of the value
    # against which to normalize. To do this, look in the xvals by experiment, 
    # then index that according to whether we cut it or not. There doesn't seem 
    # to be a cleaner way to do the indexing of cut.
    normidx = Dict("exobase"=>findfirst(isequal(meanTe), tvals["exobase"][1:1:length(tvals["exobase"])]), 
                   "tropopause"=>findfirst(isequal(meanTt), tvals["tropopause"][1:1:length(tvals["tropopause"])]), 
                   "surface"=>findfirst(isequal(meanTs), tvals["surface"][1:1:length(tvals["surface"])]),
                   )
    nomdict = makenomdict(abs_or_mr) # get info for nom case which has values not % by 10
    
    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    # Make the plot 
    fig, ax = subplots(1, 3, sharex=false, sharey=false, figsize=(21, 5))
    subplots_adjust(wspace=0.3)
    

    for i in 1:3
        plot_bg(ax[i])
        ex = exps[i]
        
        Hdiff = normalize_val(output_dict[ex]["H"], nomdict["H"])
        Ddiff = normalize_val(output_dict[ex]["D"], nomdict["D"])
        H2diff = normalize_val(output_dict[ex]["H2"], nomdict["H2"])
        HDdiff = normalize_val(output_dict[ex]["HD"], nomdict["HD"])       
        Hfdiff = normalize_val(output_dict[ex]["Hflux"], nomdict["Hflux"])
        Dfdiff = normalize_val(output_dict[ex]["Dflux"], nomdict["Dflux"])

        # plot the actual stuff
        ax[i].plot(tvals[ex], Hdiff, marker="x", ms=sz, color=c[1], zorder=10, label="H")
        ax[i].plot(tvals[ex], Ddiff, marker="*", ms=sz, color=c[2], zorder=10, label="D")
        ax[i].plot(tvals[ex], H2diff, marker="o", ms=sz, color=c[3], zorder=10, label="H2")
        ax[i].plot(tvals[ex], HDdiff, marker="D", ms=sz, color=c[4], zorder=10, label="HD")
        ax[i].plot(tvals[ex], Hfdiff, marker="^", ms=sz, color=c[5], zorder=10, label=L"\phi_H")
        ax[i].plot(tvals[ex], Dfdiff, marker="v", ms=sz, color=c[6], zorder=10, label=L"\phi_D")
        # nominal value
        ax[i].axvline(nom_i_py[ex], color=medgray, zorder=5)
        ax[i].text(gmean_txt_loc[ex][1], gmean_txt_loc[ex][2], "Global\nmean", 
                   color=medgray, ha="left", va="top")
        ax[i].axhline(1, color=medgray, zorder=5)

        # text on plots
        dfentry = linelbls[linelbls.Exp.==ex, :]
        ax[i].text(dfentry.H[1][1], dfentry.H[1][2], "H", color=c[1], ha="left", va="top")
        ax[i].text(dfentry.D[1][1], dfentry.D[1][2], "D", color=c[2], ha="left", va="top")
        ax[i].text(dfentry.H2[1][1], dfentry.H2[1][2], "H2", color=c[3], ha="left", va="top")
        ax[i].text(dfentry.HD[1][1], dfentry.HD[1][2], "HD", color=c[4], ha="left", va="top")
        ax[i].text(dfentry.fluxH[1][1], dfentry.fluxH[1][2], L"\phi_H", color=c[5], ha="left", va="top")
        ax[i].text(dfentry.fluxD[1][1], dfentry.fluxD[1][2], L"\phi_D", color=c[6], ha="left", va="top")


        # Various configurations and such 
        if ex=="surface"
            ax[i].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(20))
            ax[i].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
        elseif ex=="tropopause"
            ax[i].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
            # ax[i].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
        elseif ex=="exobase"
            ax[i].set_yscale("log")
        end

        ax[i].set_xlabel(ex*" temperature (K)")
    end
    ax[1].set_ylabel(L"X/X(\overline{T}_{global})")
    ax[1].text(150, 1.4, "a", color="black", weight="bold", va="top", 
               ha="center", fontsize=26)#
    ax[2].text(100, 1.9, "b", color="black", weight="bold", va="top", 
               ha="center", fontsize=26)
    ax[3].text(150, 100, "c", color="black", weight="bold", va="top", 
               ha="center", fontsize=26)

    # Other species that can be plotted, or species at other locations
    # HDtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ODtopdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # HDO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # DO2topdiff = normalize(MRdictvar["HDtop"], normidx[ex])
    # ax.plot(xvals[ex], HDtopdiff, marker="o", color="#7D3403", zorder=10)
    # ax.plot(xvals[ex], ODtopdiff, marker="o", color="red", zorder=10)
    # ax.plot(xvals[ex], HDO2topdiff, marker="o", color="blue", zorder=10)
    # ax.plot(xvals[ex], DO2topdiff, marker="o", color="purple", zorder=10)
    
    savefig(mainpath*"compare_nominal_temps_"*abs_or_mr*".png", bbox_inches="tight")
end

function make_T_f_plot(output_MR)
    #=
    makes a 3-panel plot of the tradeoff of f. 
    =#

    maincol = "#574D9D"
    secondcol = "#588C54"
    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    # Dictionaries for where to put things for each experiment 
    exps = ["surface", "tropopause", "exobase"]
    gmean_txt_loc = Dict("surface"=>[218, 1], "tropopause"=>[132, 1], "exobase"=>[208, 1])
    xt_args = Dict("surface"=>[150, 20, 275], "tropopause"=>[70, 30, 165], "exobase"=>[150, 50, 350])
    nom_i_py = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs)
    ylims = Dict("tropopause"=>[1e-5, 1e0], "exobase"=>[1e-5, 1e0], "surface"=>[1e-5, 1e0])
    mmin = Dict("exobase"=>25, "tropopause"=>10, "surface"=>10)
    mmaj = Dict("exobase"=>50, "tropopause"=>30, "surface"=>20)

    # Make the figure 
    fig, ax = subplots(1, 3, sharex=false, sharey=false, figsize=(21, 5))
    subplots_adjust(wspace=0.3)
    
    for i in 1:3
        plot_bg(ax[i])
        ex = exps[i]

        # plot
        ax[i].plot(tvals[ex], output_MR[ex]["f"], marker="o", color=maincol, zorder=10)
        ax[i].axvline(nom_i_py[ex], color=medgray, zorder=5)
        ax[i].text(gmean_txt_loc[ex][1], gmean_txt_loc[ex][2], "Global\nmean", 
                   color=medgray, ha="left", va="top")

        # axis format
        mpltic = pyimport("matplotlib.ticker")
        ax[i].set_xlabel(ex*" temperature (K)", fontsize=24)
        ax[i].xaxis.set_major_locator(mpltic.MultipleLocator(mmaj[ex]))
        ax[i].xaxis.set_minor_locator(mpltic.MultipleLocator(mmin[ex]))
        ax[i].set_xticks(xt_args[ex][1]:xt_args[ex][2]:xt_args[ex][3], minor=false)
        ax[i].tick_params(axis="y", labelcolor=maincol, which="both")
        ax[i].tick_params(axis="x", which="both")
        # setp(ax.get_xticklabels(), fontsize=22)

        ax[i].set_yscale("log") 
        ax[i].set_ylim(ylims[ex][1], ylims[ex][2])
        ax[i].set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
        # setp(ax.get_yticklabels(), fontsize=22)

        # plot the secondary axis showing increase
        nomdict = makenomdict("mr")
        fdiff = normalize_val(output_MR[ex]["f"], nomdict["f"])
        ax_2 = ax[i].twinx()
        if i==3
            ax_2.set_ylabel(L"$f$/$f(\overline{T}_{global})$", color=secondcol, fontsize=24)
        end
        ax_2.plot(tvals[ex], fdiff, marker="o", color=secondcol)
        ax_2.tick_params(axis="y", labelcolor=secondcol, which="both")
        for side in ["top", "bottom", "left", "right"]
            ax_2.spines[side].set_visible(false)
        end
    end


    
    ax[1].set_ylabel(L"Fractionation factor $f$", color=maincol, fontsize=24) 
    ax[1].text(150, 1, "a", color="black", fontsize=26, weight="bold", va="top")
    ax[2].text(100, 1, "b", color="black", fontsize=26, weight="bold", va="top")
    ax[3].text(150, 1, "c", color="black", fontsize=26, weight="bold", va="top")


    savefig(mainpath*"f_vs_temps.png", bbox_inches="tight")
end

function make_T_output_vs_data(output_MR, output_abs)
    #=
    Makes the plots comparing simulation output with observational data.

    =#

    # set up
    exps = ["surface", "tropopause", "exobase"]
    c1 = ["#10007A", "#2F7965", "#e16262", "#D59A07"]  # just colors
    medgray = "#444444"

    nom_i_py = Dict("exobase"=>meanTe, "tropopause"=>meanTt, "surface"=>meanTs)


    # where to place text, basically.
    mean_text = Dict("surface"=>[218, 100], "tropopause"=>[131, 4], "exobase"=>[208, 4])
    ylims = Dict("surface"=>[-10, 190], "tropopause"=>[-6.3, 4.2], "exobase"=>[-6.3, 4.2])
    xt_args = Dict("surface"=>150:20:280, "tropopause"=>100:10:160, "exobase"=>150:50:350)
    linelbls = DataFrame(Exp=["surface", "tropopause", "exobase"], 
                          CO=[[155, -3],  [100, -5],   [175, -5]], 
                          O2=[[145, 40],  [115, -2],   [150, -2]], 
                          H2=[[180, 1],   [100, 2.3],  [150, 4.1]],
                          O3=[[155, 110], [100, -0.8], [155, -0.8]])

    # make plots pretty
    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 22
    rcParams["axes.labelsize"]= 24
    rcParams["xtick.labelsize"] = 22
    rcParams["ytick.labelsize"] = 22

    fig, ax = subplots(1, 3, sharex=false, sharey=false, figsize=(21, 5))
    subplots_adjust(wspace=0.15)


    for i in 1:3
        plot_bg(ax[i])
        ex = exps[i]

        # calculate the relativeness of each point wrt data
        COdiff = (output_MR[ex]["CO"] .- data[1])/s[1]
        O2diff = (output_MR[ex]["O2"] .- data[2])/s[2]
        # H2diff = (ABSdictvar["H2"] .- data[3])/s[3]  # absolute abundance
        H2diff = (output_MR[ex]["H2MR"]./1e-6 .- data[3])/s[3] # ppm
        O3diff = (areadensity_to_micron_atm(output_abs[ex]["O3"]) .- data[4])/s[4] 

        # plot
        ax[i].plot(tvals[ex], COdiff, marker="o", ms=sz, color=c1[1], zorder=10)
        ax[i].plot(tvals[ex], O2diff, marker="x", ms=sz, color=c1[2], zorder=10)
        ax[i].plot(tvals[ex], H2diff, marker="*", ms=sz, color=c1[3], zorder=10)
        ax[i].plot(tvals[ex], O3diff, marker="^", ms=sz, color=c1[4], zorder=10)
        # plot the mean values
        ax[i].axvline(nom_i_py[ex], color=medgray, zorder=5)
        ax[i].axhline(0, color="black")

        # text
        ax[i].text(mean_text[ex][1], mean_text[ex][2], "Global\nmean", color=medgray, ha="left", va="top")
        dfentry = linelbls[linelbls.Exp.==ex, :]
        ax[i].text(dfentry.CO[1][1], dfentry.CO[1][2], "CO", color=c1[1], ha="left", va="top")
        ax[i].text(dfentry.O2[1][1], dfentry.O2[1][2], L"O$_2$", color=c1[2], ha="left", va="top")
        ax[i].text(dfentry.H2[1][1], dfentry.H2[1][2], L"H$_2$", color=c1[3], ha="left", va="top")
        ax[i].text(dfentry.O3[1][1], dfentry.O3[1][2], L"O$_3$", color=c1[4], ha="left", va="top")

        # labels and such
        ax[i].set_xlabel(ex*" temperature (K)")
        ax[i].set_xticks(xt_args[ex])
        ax[i].set_ylim(ylims[ex][1], ylims[ex][2])
    end
   

    # SURFACE PLOT CONFIG
    ax[1].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    ax[1].set_yscale("symlog", linthresh=0.2)
    ax[1].set_yticks([-10, -1, 0, 1, 10, 100])
    ax[1].set_yticklabels([-10, -1, 0, 1, 10, 100])
    ax[1].set_ylabel(L"($X_{model}-X_{obs}$)/$\sigma$")
    
    # TROPOPAUSE PLOT CONFIG
    ax[2].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    ax[2].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
    
    # Extra clarifying text
    ax[3].text(365, 4, "Model value \n"*L"$>$ observations", va="top")
    ax[3].text(365, 0, "Model value \n"*L"$\approx$observations", va="center")
    ax[3].text(365, -4, "Model value \n"*L"$<$observations")

    ax[1].text(270, 80, "a", color="black", fontsize=26, weight="bold", va="bottom")
    ax[2].text(160, 4, "b", color="black", fontsize=26, weight="bold", va="top")
    ax[3].text(350, 4, "c", color="black", fontsize=26, weight="bold", va="top")
    
    savefig(mainpath*"output_vs_data_temps.png", bbox_inches="tight")      
end

# Setup and folder Location ====================================================
# mainpath = "/data/GDrive-CU/Research/Results/TradeoffPlots/Tradeoffs - solar mean/"
mainpath = "/home/emc/GDrive-CU/Research/Results/"
println("Enter a folder to use, no slashes (/home/emc/GDrive-CU/Research/Results/<folder>: ")
append_me = readline(stdin)
mainpath = mainpath*append_me*"/"

makeplots = false   # whether to make the plots that go with each experiment
other_deuterated = false
write_new_files = false  # set to true if running for first time after new simulations

# Function calls ===============================================================

println("Analyzing water model output, building dicts, making supporting plots")
water_data_mr = analyze_water("mr", other_deuterated, makeplots)
water_data_abs = analyze_water("abs", other_deuterated, makeplots)

# println()
# println("Analyzing O flux model output, building dicts, making supporting plots")
# o_data_abs = analyze_Oflux("abs", other_deuterated, makeplots)
# o_data_mr = analyze_Oflux("mr", other_deuterated, makeplots)

# println()
println("Analyzing temp model output, building dicts, making supporting plots")
T_data_mr = analyze_T("mr", other_deuterated, makeplots)
T_data_abs = analyze_T("abs", other_deuterated, makeplots)

println("Writing jld storage files")
if write_new_files
    wd_mr = jldopen(mainpath*"water_MR_data.jld", "w")
    @write wd_mr water_data_mr
    close(wd_mr)
    wd_abs = jldopen(mainpath*"water_abs_data.jld", "w")
    @write wd_abs water_data_abs
    close(wd_abs)

    # O_mr = jldopen(result_folder*"O_MR_data.jld", "w")
    # @write O_mr o_data_mr
    # close(O_mr)
    # O_abs = jldopen(result_folder*"O_abs_data.jld", "w")
    # @write O_abs o_data_abs
    # close(O_abs)

    T_mr = jldopen(mainpath*"T_MR_data.jld", "w")
    @write T_mr T_data_mr
    close(T_mr)
    T_abs = jldopen(mainpath*"T_abs_data.jld", "w")
    @write T_abs T_data_abs
    close(T_abs)
end

# Result and discussion plots: =================================================
# 1. f as a function of parameter (e.g. temperature or water vapor)
# 2. atomic/molecular H/D and fluxes compared across the simulations
# 3. Comparison of model output of CO, O2, O3, H2 to measurements in lit

# make_rel_change_plots(o_data_mr, o_data_abs, "O flux")

# Temperature plots
println("Making temperature plots")
# make_T_f_plot(T_data_mr)
# make_T_Hspecies_plot(T_data_mr, "mr")
# make_T_Hspecies_plot(T_data_abs, "abs")
make_T_output_vs_data(T_data_mr, T_data_abs)

# Water plots
# println("Making water plots")
# make_water_f_plot(water_data_mr)
# make_water_Hspecies_plot(water_data_mr, "mr")
# make_water_Hspecies_plot(water_data_abs, "abs")
# make_water_output_vs_data(water_data_mr, water_data_abs)