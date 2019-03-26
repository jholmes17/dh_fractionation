################################################################################
# plot_snapshot.jl - Plots a snapshot of the atmospheric state at a specified
# time index within the simulation.
#
# Eryn Cangi
# 31 August 2018
# Currently tested for Julia: 0.7
################################################################################

using PyPlot
using HDF5
using LaTeXStrings

const fullspecieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
                         :O1D, :H, :N2, :Ar, :CO2pl, :HOCO,
                         # species for deuterium chemistry:
                         :HDO, :OD, :HDO2, :D, :DO2, :HD, :DOCO];
specieslist=fullspecieslist;

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

function n_tot(n_current, z)
    # get the total number density at a given altitude
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

function plot_snap(readfile, timevec, poster, convergefile=false, init=false)
    #=
    Creates a plot of selected species volume mixing ratio by altitude at the
    equilibrium case. The data is given in number density by species and
    altitude, so this also converts number density to volume mixing ratio first.
    Plots only the species:

    readfile: hdf5 file of Results
    timevec: time to plot, must be a 1-length vector so I dont have to change code lol
    poster: make fonts big if this is going on a poster
    convergefile: whether the file being read is a convergence (equilibrium) file
    init: whether the snapshot wanted is the beginning of the simulation. init
          must be false if convergefile is true

    =#
    println(readfile)
    const alt = h5read(readfile,"n_current/alt")  # read in the altitudes
    if convergefile != true
        const times = h5read(readfile,"n_current/timelist")  # read in the times
    end

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

    # Font size settings
    if poster==true
        fs = Dict("ticks"=>16, "labels"=>22, "legend"=>18, "title"=>26, "spec"=>14)
    elseif poster==false
        fs = Dict("ticks"=>16, "labels"=>20, "legend"=>14, "title"=>22, "spec"=>14)
    end

    # make a figure and set up some global labels
    fig, ax = subplots(figsize=(10,6))

    # do not change the following line and do not use tight_layout(), it conflicts
    subplots_adjust(wspace=0, bottom=0.15)
    ax[:set_xlabel]("Volume Mixing Ratio [ppm]", fontsize=fs["labels"])
    ax[:set_ylabel]("Altitude [km]", fontsize=fs["labels"])
    y = alt[2:end-1]./10^5  # altitudes are always the same

    # Loop over times to make subplots
    for i in (1:length(timevec))
        # grab the data from just the iteration we want
        if init==false
            if convergefile == true
                ncur_mat = h5read(readfile, "n_current/n_current_mat")
            else
                ncur_mat = h5read(readfile,"n_current/iter_$(timevec[i])");
            end
        elseif init==true
            ncur_mat = h5read(readfile,"n_current/init");
        end

        # reformat the matrix into an easy-to-read dictionary
        ncurrent = Dict{Symbol,Array{Float64, 1}}()
        for s in [1:length(species);]  # the ; makes it not concatenate
            ncurrent[species[s]] = reshape(ncur_mat[:,s],length(alt)-2)
        end

        # convert to volume mixing ratio
        ncurrent_vmr = deepcopy(ncurrent)

        for k in keys(ncurrent_vmr)
            ncurrent_vmr[k] = ncurrent_vmr[k] ./ map(z->n_tot(ncurrent, z), alt[2:end-1])
        end

        # convert the iteration number to a time in years, days, hours
        timestr = ""
        if convergefile == false
            secs = times[timevec[i]]
            time = convert_secs(secs)
            key = ["years ", "days ", "hours "]

            for t in range(1,length(time))
              if time[t] != 0
                  timestr = timestr * "$(Int(time[t])) $(key[t])"
              end
            end
            if timestr == ""
                timestr = "Initial"
            end
        else
            timestr = "Converged atomsphere"
        end

        # PLOT - SELECT SPECIES
        for (s, sstr) in zip(species, species_str)
            # only plot some of the most important species
            if !contains(string(s), "J")  # do not plot J rates
                if sstr in ["O_3", "H_2O_2", "HDO_2", "DO_2", "HO_2", "HD", "H_2", "H", "D", "H_2O", "HDO", "O", "O_2", "CO_2"]
                    x = ncurrent_vmr[s]

                    # set line thicknesses
                    if sstr in ["H", "D"]
                        lw = 3
                    else
                        lw = 1
                    end
                    ax[:semilogx](x, y, linewidth=lw, alpha=0.85,
                                       color=speciescolor[s],
                                       linestyle=speciesstyle[s])

                    # stagger labels between top and bottom of plot =, easier to read
                    if sstr in ["HDO_2", "H_2O_2", "D", "H_2", "O_2"]
                        ax[:text](0.4*x[1], y[1]-12, latexstring(sstr),
                                       color=speciescolor[s], fontsize=fs["spec"])
                    elseif sstr=="CO_2"
                        ax[:text](0.1*x[1], y[1]-12, latexstring(sstr),
                                       color=speciescolor[s], fontsize=fs["spec"])
                    else
                        ax[:text](0.4*x[end], y[end]+4, latexstring(sstr),
                                       color=speciescolor[s], fontsize=fs["spec"])
                    end

                end

                # set the xlim, ylim
                ax[:set_ylim](-15,215)
                ax[:set_xticks]([1e0, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15])
                #ax[:set_xlim](1e-17,5)
                ax[:tick_params]("both",labelsize=fs["ticks"])

            end
        end
    end # loop over time vector
    savefig("../notebookpics/snapshot_MAYBEFIX"*string(timevec[1])*".png")
end


# base = "/data/GDrive-CU/Research/Results/Mikes_Results"
# rf = base*"/one_year_response_to_80ppm_at_60km.h5"
# plot_snap(rf, [1], true, false, false)


base = "/data/GDrive-CU/Research/Results/VarWaterTemp"
rf = base*"/temp_192_110_199/converged_standardwater_D_temp_192_110_199.h5"
const alt = h5read(rf,"n_current/alt")
n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])


plot_snap(rf, [1299], true, true, false) # for convergence files
