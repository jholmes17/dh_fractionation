using HDF5
using PyPlot
using LaTeXStrings

function convert_secs(secs)
    #=
    A function to convert a given amount of seconds into a sensibly rounded
    number of years, days, and hours. Since the simulation only goes for a year,
    we don't care if there are extra days or hours past one year, and so we just
    truncate.
    =#

    # using the div function gets us the quotient as an integer.
    yrs = div(secs, 3.154E7)
    days = div((secs % 3.154E7), 86400)
    hrs = div((secs % (3.154E7)) % 86400, 3600)
    mins = div(((secs % (3.154E7)) % 86400) % 3600, 60)

    # toggle these print statements on and off if you need them
    # println("Years: $(yrs)")
    # println("Days: $(days)")
    # println("hours: $(hrs)")
    # println("mins: $(mins)")

    # these limits are sort of arbitrary. 12+ hours is enough to round up to the
    # next day, but i don't believe that 180 days is enough to round to another
    # year.
    if mins > 30
        hrs += 1
        mins = 0
    end
    if hrs > 12
        days += 1
        hrs = 0
    end
    if days > 320
        yrs += 1
        days = 0
    end
    if yrs >= 1  # at this time scale, we don't care how many days or hours
        days = 0
        hrs = 0
    end
    return [yrs, days, hrs]
end

function plot_select_spec(readfile, timevec, poster, ext="")
    #=
    Creates a plot of selected species volume mixing ratio by altitude at the
    equilibrium case. The data is given in number density by species and
    altitude, so this also converts number density to volume mixing ratio first.
    Plots only the species:

    =#

    const alt = h5read(readfile,"n_current/alt")  # read in the altitudes
    const times = h5read(readfile,"n_current/timelist")  # read in the times

    # read in species names and map them to symbols
    species = map(Symbol,h5read(readfile,"n_current/species"))
    species_str = map(String,h5read(readfile,"n_current/species"))
    # make some latex-formatted labels
    for n in (1:length(species))
        if species_str[n][1] != 'J'
            species_str[n] = replace(species_str[n], '3', "_3")
            species_str[n] = replace(species_str[n], '2', "_2")
            species_str[n] = replace(species_str[n], "1D", "(^1D)")
            species_str[n] = replace(species_str[n], "pl", "^+")
        end
    end

    # colors
    speciescolor = Dict(
        # H group
        :H => "#ff0000", # red
        :D => "#ff0000", #
        :H2 => "#e526d7", #
        :HD =>  "#e526d7", # dark pink/magenta

        # hydroxides
        :OH => "#7700d5",
        :OD => "#7700d5", # purple

        # water group (roughly, I ain't a chemist)
        :H2O => "#0083dc",  # cornflower blue
        :HDO => "#0083dc",
        :H2O2 => "#0000ff",
        :HDO2 => "#0000ff", # true blue
        :HO2 => "#046868",
        :DO2 => "#046868",  # dark teal

        # O group
        :O1D => "#808000",  # olive
        :O => "#1a6115",    # forest green
        :O2 => "#15da09",    # kelly/grass green
        :O3 => "#269e56",   # light green

        # CO group
        :CO2 => "#d18564",   # dark peach
        :CO2pl => "#614215", # brown
        :CO => "#ff6600",    # orange
        :HOCO => "#e8ba8c", #tannish
        :DOCO => "#e8ba8c",

        # nonreactants
        :Ar => "#808080",
        :N2 => "#cccccc",);

    speciesstyle = Dict(
        # H group
        :H => "-",
        :D => "--",
        :H2 => "-",
        :HD => "--",

        # hydroxides
        :OH => "-", # Lilac
        :OD => "--", # purple

        # water group (roughly, I ain't a chemist)
        :H2O => "-",  # cornflower blue
        :HDO => "--",  # navy
        :H2O2 => "-", # nice blue
        :HDO2 => "--", # true blue
        :HO2 => "-",  # teal
        :DO2 => "--",  # dark teal

        # O group
        :O1D => "-",  # olive
        :O => "-",    # forest green
        :O2 => "-",    # kelly/grass green
        :O3 => "-",   # light green

        # CO group
        :CO2 => "-",   # dark peach
        :CO2pl => "-", # brown
        :CO => "-",    # orange
        :HOCO => "-", #tannish
        :DOCO => "--",

        # nonreactants
        :Ar => "-",
        :N2 => "-",);

    # Font size setup_photochemistry
    if poster==true
        fs = Dict("ticks"=>16, "labels"=>22, "legend"=>18, "title"=>26, "spec"=>14)
    elseif poster==false
        fs = Dict("ticks"=>16, "labels"=>20, "legend"=>14, "title"=>22, "spec"=>14)
    end

    # make a figure and set up some global labels
    fig, ax = subplots(1,length(timevec), sharey=true, figsize=(25,5))

    # do not change the following and do not use tight_layout(), it conflicts 
    subplots_adjust(wspace=0, bottom=0.15)
    ax[1,1][:set_xlabel]("Volume Mixing Ratio [ppm]", fontsize=fs["labels"])
    ax[1,1][:set_ylabel]("Altitude [km]", fontsize=fs["labels"])
    y = alt[2:end-1]./10^5  # altitudes are always the same
    # stuff to put on return plots
    returntext = ["", "", "", ""]

    # Loop over times to make subplots
    for i in (1:length(timevec))
        # grab the data from just the iteration we want
        ncur_mat = h5read(readfile,"n_current/iter_$(timevec[i])");

        # reformat the matrix into an easy-to-read dictionary
        ncurrent = Dict{Symbol,Array{Float64, 1}}()
        for s in [1:length(species);]  # the ; makes it not concatenate
            ncurrent[species[s]] = reshape(ncur_mat[:,s],length(alt)-2)
        end

        # convert to volume mixing ratio
        ncurrent_vmr = deepcopy(ncurrent)

        # iterate over altitudes
        for a in range(1,length(alt)-2)
            n_tot_this_alt = 0

            # including the following 2 loops separately is less code than
            # copying in n_tot and all its dependencies from setup_photochemistry

            # compute the total number density at altitude alt[a]
            for s in species
                n_tot_this_alt += ncurrent[s][a]
             end

            # convert to VMR and fill the array
            for s2 in species
                ncurrent_vmr[s2][a] = ncurrent[s2][a] / n_tot_this_alt
            end
        end

        # convert the iteration number to a time in years, days, hours
        secs = times[timevec[i]]
        time = convert_secs(secs)
        key = ["years ", "days ", "hours "]
        timestr = ""
        for t in range(1,length(time))
          if time[t] != 0
              timestr = timestr * "$(Int(time[t])) $(key[t])"
          end
        end
        # handle the case where the simulation is in the initial state (time = 0 approx)
        if timestr == ""
            timestr = "Initial"
        end

        # PLOT - SELECT SPECIES
        for (s, sstr) in zip(species, species_str)
            # only plot some of the most important species
            if !contains(string(s), "J")  # do not plot J rates
                if sstr in ["H_2O_2", "HDO_2", "DO_2", "HO_2", "HD", "H_2", "H", "D", "H_2O", "HDO", "O", "O_2", "CO_2"]
                    x = ncurrent_vmr[s]
                    # set line thicknesses
                    if sstr in ["H", "D"]
                        lw = 3
                    else
                        lw = 1
                    end
                    ax[i,1][:semilogx](x, y, linewidth=lw, alpha=0.85,
                                       color=speciescolor[s],
                                       linestyle=speciesstyle[s])

                    # stagger labels between top and bottom of plot =, easier to read
                    if sstr in ["HDO_2", "H_2O_2", "D", "H_2", "O_2"]
                        ax[i,1][:text](0.4*x[1], y[1]-12, latexstring(sstr),
                                       color=speciescolor[s], fontsize=fs["spec"])
                    elseif sstr=="CO_2"
                        ax[i,1][:text](0.1*x[1], y[1]-12, latexstring(sstr),
                                       color=speciescolor[s], fontsize=fs["spec"])
                    else
                        ax[i,1][:text](0.4*x[end], y[end]+4, latexstring(sstr),
                                       color=speciescolor[s], fontsize=fs["spec"])
                    end

                end

                # title stuff
                if ext==""
                    ax[i,1][:set_title](timestr, fontsize=fs["title"])
                else
                    ax[i,1][:text](10.0^-17, 180, returntext[i])
                end

                # set the xlim, ylim
                ax[i,1][:set_ylim](-15,215)
                ax[i,1][:set_xticks]([1e0, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15])
                #ax[i,1][:set_xlim](1e-17,5)
                ax[i,1][:tick_params]("both",labelsize=fs["ticks"])

            end
        end
    end # loop over time vector
    savefig("1yr_evolution"*ext*".png")
end

#TODO: fix this so that it is like select_spec to ensure best plotting
function plot_all_spec(readfile, timevec, ext="")
    #=
    Creates a plot of all species volume mixing ratio by altitude (altitude is
    on the y axis as is traditional) at the equilibrium case. The data is given
    in number density by species and altitude, so this also converts number
    density to volume mixing ratio first.
    =#

    const alt = h5read(readfile,"n_current/alt")  # read in the altitudes
    const times = h5read(readfile,"n_current/timelist")  # read in the times

    # read in species names and map them to symbols
    species = map(Symbol,h5read(readfile,"n_current/species"))
    species_str = map(String,h5read(readfile,"n_current/species"))
    # make some latex-formatted labels
    for n in (1:length(species))
        if species_str[n][1] != 'J'
            species_str[n] = replace(species_str[n], '3', "_3")
            species_str[n] = replace(species_str[n], '2', "_2")
            species_str[n] = replace(species_str[n], "1D", "(^1D)")
            species_str[n] = replace(species_str[n], "pl", "^+")
        end
    end

    # colors
    speciescolor = Dict(
        # H group
        :H => "#ff0000", # red
        :D => "#ff0000", #
        :H2 => "#e526d7", #
        :HD =>  "#e526d7", # dark pink/magenta

        # hydroxides
        :OH => "#7700d5",
        :OD => "#7700d5", # purple

        # water group (roughly, I ain't a chemist)
        :H2O => "#0083dc",  # cornflower blue
        :HDO => "#0083dc",
        :H2O2 => "#0000ff",
        :HDO2 => "#0000ff", # true blue
        :HO2 => "#046868",
        :DO2 => "#046868",  # dark teal

        # O group
        :O1D => "#808000",  # olive
        :O => "#1a6115",    # forest green
        :O2 => "#15da09",    # kelly/grass green
        :O3 => "#269e56",   # light green

        # CO group
        :CO2 => "#d18564",   # dark peach
        :CO2pl => "#614215", # brown
        :CO => "#ff6600",    # orange
        :HOCO => "#e8ba8c", #tannish
        :DOCO => "#e8ba8c",

        # nonreactants
        :Ar => "#808080",
        :N2 => "#cccccc",);
    speciesstyle = Dict(
        # H group
        :H => "-",
        :D => "--",
        :H2 => "-",
        :HD => "--",

        # hydroxides
        :OH => "-", # Lilac
        :OD => "--", # purple

        # water group (roughly, I ain't a chemist)
        :H2O => "-",  # cornflower blue
        :HDO => "--",  # navy
        :H2O2 => "-", # nice blue
        :HDO2 => "--", # true blue
        :HO2 => "-",  # teal
        :DO2 => "--",  # dark teal

        # O group
        :O1D => "-",  # olive
        :O => "-",    # forest green
        :O2 => "-",    # kelly/grass green
        :O3 => "-",   # light green

        # CO group
        :CO2 => "-",   # dark peach
        :CO2pl => "-", # brown
        :CO => "-",    # orange
        :HOCO => "-", #tannish
        :DOCO => "--",

        # nonreactants
        :Ar => "-",
        :N2 => "-",);
    # make a figure and set up some global labels
    fig, ax = subplots(1,4, sharey=true, figsize=(18,5))
    subplots_adjust(wspace=0)
    ax[1,1][:set_xlabel]("Volume Mixing Ratio [ppm]")
    ax[1,1][:set_ylabel]("Altitude [km]")
    y = alt[2:end-1]./10^5  # altitudes are always the same
    # stuff to put on return plots
    returntext = ["", "", "", ""]

    # Loop over times to make subplots
    for i in (1:length(timevec))
        # grab the data from just the iteration we want
        ncur_mat = h5read(readfile,"n_current/iter_$(timevec[i])");

        # reformat the matrix into an easy-to-read dictionary
        ncurrent = Dict{Symbol,Array{Float64, 1}}()
        for s in [1:length(species);]  # the ; makes it not concatenate
          ncurrent[species[s]] = reshape(ncur_mat[:,s],length(alt)-2)
        end

        # convert to volume mixing ratio
        ncurrent_vmr = deepcopy(ncurrent)

        # iterate over altitudes
        for a in range(1,length(alt)-2)
          n_tot_this_alt = 0

          # including the following 2 loops separately is less code than
          # copying in n_tot and all its dependencies from setup_photochemistry

          # compute the total number density at altitude alt[a]
          for s in species
              n_tot_this_alt += ncurrent[s][a]
          end

          # convert to VMR and fill the array
          for s2 in species
              ncurrent_vmr[s2][a] = ncurrent[s2][a] / n_tot_this_alt
          end
        end

        # convert the iteration number to a time in years, days, hours
        secs = times[timevec[i]]
        time = convert_secs(secs)
        key = ["years ", "days ", "hours "]
        timestr = ""
        for t in range(1,length(time))
          if time[t] != 0
              timestr = timestr * "$(Int(time[t])) $(key[t])"
          end
        end
        # handle the case where the simulation is in the initial state (time = 0 approx)
        if timestr == ""
          timestr = "Initial"
        end

        # PLOT - ALL SPECIES
        for (s, sstr) in zip(species, species_str)
            # only plot some of the most important species
            if !contains(string(s), "J")  # do not plot J rates
                x = ncurrent_vmr[s]

                if sstr in ["HDO_2", "H_2O_2", "DO_2", "HO_2", "H", "D"]
                    lw = 3
                else
                    lw = 1
                end

                ax[i,1][:semilogx](x, y, linewidth=1, alpha=0.85,
                                   color=speciescolor[s],
                                   linestyle=speciesstyle[s], linewidth=lw)
                # stagger labels between top and bottom of plot =, easier to read
                if sstr in ["HDO_2", "H_2O_2", "HD", "H_2", "H", "CO_2", "CO", "Ar", "O(^1D)", "O"]
                    ax[i,1][:text](0.2*x[1], y[1]-9, latexstring(sstr),
                                   color=speciescolor[s])
                else
                    ax[i,1][:text](0.2*x[end], y[end]+4, latexstring(sstr),
                                   color=speciescolor[s])
                end

                # title stuff
                if ext==""
                    ax[i,1][:set_title](timestr)
                else
                    ax[i,1][:text](10.0^-17, 180, returntext[i])
                end
            end
        end
    end # loop over time vector
    ylim(-10,210)
    tight_layout()
    savefig("1yr_evolution"*ext*".png", bbox_to_anchor="tight")
end

#TODO: fix this so that it is like select_spec to ensure best plotting
function plot_top_spec(readfile, timevec, ext="")
    #=
    Creates a plot of the most populous species volume mixing ratio by altitude
    at the equilibrium case.

    timevec: vector of times at which to make plots
    =#
    const alt = h5read(readfile,"n_current/alt")  # read in the altitudes
    const times = h5read(readfile,"n_current/timelist")  # read in the times

    # read in species names and map them to symbols
    species = map(Symbol,h5read(readfile,"n_current/species"))
    species_str = map(String,h5read(readfile,"n_current/species"))
    # make some latex-formatted labels
    for n in (1:length(species))
        if species_str[n][1] != 'J'
            species_str[n] = replace(species_str[n], '3', "_3")
            species_str[n] = replace(species_str[n], '2', "_2")
            species_str[n] = replace(species_str[n], "1D", "(^1D)")
            species_str[n] = replace(species_str[n], "pl", "^+")
        end
    end

    # colors: same dict as in setup_photochemistry.
    speciescolor = Dict(
        # H group
        :H => "#ff0000", # red
        :D => "#ff0000", #
        :H2 => "#e526d7", #
        :HD =>  "#e526d7", # dark pink/magenta

        # hydroxides
        :OH => "#7700d5",
        :OD => "#7700d5", # purple

        # water group (roughly, I ain't a chemist)
        :H2O => "#0083dc",  # cornflower blue
        :HDO => "#0083dc",
        :H2O2 => "#0000ff",
        :HDO2 => "#0000ff", # true blue
        :HO2 => "#046868",
        :DO2 => "#046868",  # dark teal

        # O group
        :O1D => "#808000",  # olive
        :O => "#1a6115",    # forest green
        :O2 => "#15da09",    # kelly/grass green
        :O3 => "#269e56",   # light green

        # CO group
        :CO2 => "#d18564",   # dark peach
        :CO2pl => "#614215", # brown
        :CO => "#ff6600",    # orange
        :HOCO => "#e8ba8c", #tannish
        :DOCO => "#e8ba8c",

        # nonreactants
        :Ar => "#808080",
        :N2 => "#cccccc",);

    # make a figure and set up some global labels
    fig, ax = subplots(1,4, sharey=true, figsize=(18,5))
    ax[1,1][:set_xlabel]("Volume Mixing Ratio [ppm]")
    ax[1,1][:set_ylabel]("Altitude [km]")
    y = alt[2:end-1]./10^5  # altitudes are always the same
    # stuff to put on return plots
    returntext = ["", "", "", ""]

    # Loop over times to make subplots
    for i in (1:length(timevec))# #
        # grab the data from just the iteration we want
        ncur_mat = h5read(readfile,"n_current/iter_$(timevec[i])");

        # reformat the matrix into an easy-to-read dictionary
        ncurrent = Dict{Symbol,Array{Float64, 1}}()
        for s in [1:length(species);]  # the ; makes it not concatenate
            ncurrent[species[s]] = reshape(ncur_mat[:,s],length(alt)-2)
        end

        # convert to volume mixing ratio
        ncurrent_vmr = deepcopy(ncurrent)

        # iterate over altitudes
        for a in range(1,length(alt)-2)
            n_tot_this_alt = 0

            # including the following 2 loops separately is less code than
            # copying in n_tot and all its dependencies from setup_photochemistry

            # compute the total number density at altitude alt[a]
            for s in species
                n_tot_this_alt += ncurrent[s][a]
            end

            # convert to VMR and fill the array
            for s2 in species
                ncurrent_vmr[s2][a] = ncurrent[s2][a] / n_tot_this_alt
            end
        end

        # convert the iteration number to a time in years, days, hours
        secs = times[timevec[i]]
        time = convert_secs(secs)
        key = ["years ", "days ", "hours "]
        timestr = ""
        for t in range(1,length(time))
            if time[t] != 0
                timestr = timestr * "$(Int(time[t])) $(key[t])"
            end
        end
        # handle the case where the simulation is in the initial state (time = 0 approx)
        if timestr == ""
            timestr = "Initial"
        end

        # Plotting
        for (s, sstr) in zip(species, species_str)
            # only plot some of the most important species
            if sstr in ["HDO", "H_2O", "H", "D", "HD", "H_2", "CO_2", "CO", "O", "O_2"]
                x = ncurrent_vmr[s]
                ax[i,1][:semilogx](x, y, linewidth=1, alpha=0.85, color=speciescolor[s])

                # place the species labels
                if sstr in ["HD", "HDO", "CO", "H_2O"]
                    ax[i,1][:text](0.2*x[end], y[end]+4, latexstring(sstr), color=speciescolor[s])
                else
                    ax[i,1][:text](0.2*x[1], y[1]-9, latexstring(sstr), color=speciescolor[s])
                end

                # title stuff
                if ext==""
                    ax[i,1][:set_title](timestr)
                else
                    ax[i,1][:text](10.0^-17, 180, returntext[i])
                end

                ax[i,1][:grid]("on")
            end
        end
    end # loop over time vector
    ylim(-10,210)
    savefig("1yr_evolution"*ext*".png")
end

lead = "/home/emc/Google Drive/"#"/data/GoogleDrive/"#
# f1 = "/home/emc/Google Drive/Phys/LASP/Mars/chaffincode-working/Results-Standard Water/one_year_response_to_80ppm_at_60km.h5"
f1 = lead*"Phys/LASP/Mars/chaffincode-working/one_year_response_to_80ppm_at_60km.h5"
plot_select_spec(f1, [1, 606, 790, 999], true)
#plot_all_spec(f1, [1, 606, 790, 999], true)

f2 = lead*"Phys/LASP/Mars/chaffincode-working/one_year_response_to_80ppm_at_60km_return.h5"
plot_select_spec(f2, [1, 606, 790, 999], true, "_return")
#plot_all_spec(f2, [1, 606, 790, 999], true, "_return")

# 10 days is 790 if you need it
