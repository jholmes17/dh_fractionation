using PyPlot
using HDF5

function conc_over_time_water(rf)
    #=
    returns arrays of concentrations of H2O, HDO at the ground (column 1) over time.

    rf: file containing concentrations over time (main simulation results)
    =#
    conc0_h2o = Array{Float64}(1299);  # 1299 is the number of timesteps. change as needed
    conc0_hdo = Array{Float64}(1299);

    for i in range(1, 1299) # this is all the times

        # identify the row numbers for H2O and HDO
        const slist = h5read(rf, "n_current/species") # get the list of species
        h2o_row_idx = findin(slist, ["H2O"])[1]
        hdo_row_idx = findin(slist, ["HDO"])[1]

        # get hte concentrations of H2O, HDO at ground level
        const concs = h5read(rf, "n_current/iter_"*string(i))
        conc0_h2o[i] = concs[1, h2o_row_idx]# 1st column => (alt = 0)
        conc0_hdo[i] = concs[1, hdo_row_idx]# 1st column => (alt = 0)
    end

    return conc0_h2o, conc0_hdo
end

function conc_over_time_molecH(rf)
    #=
    returns arrays of concentrations of H2, HD at the ground (column 1) over time.

    rf: file containing concentrations over time (main simulation results)
    =#
    conc0_h2 = Array{Float64}(1299);  # 1299 is the number of timesteps. change as needed
    conc0_hd = Array{Float64}(1299);

    for i in range(1, 1299) # this is all the times

        # identify the row numbers for H2O and HDO
        const slist = h5read(rf, "n_current/species") # get the list of species
        h2_row_idx = findin(slist, ["H2"])[1]
        hd_row_idx = findin(slist, ["HD"])[1]

        # get hte concentrations of H2O, HDO at ground level
        const concs = h5read(rf, "n_current/iter_"*string(i))
        conc0_h2[i] = concs[1, h2_row_idx]# 1st column => (alt = 0)
        conc0_hd[i] = concs[1, hd_row_idx]# 1st column => (alt = 0)
    end

    return conc0_h2, conc0_hd
end

function calcFractionation(rfH, rfD, p, titleext="", poster=false)
    #=
    Creates a plot of R over time for each ppm and altitude.

    rfH: readfile for H escape flux
    rfD: readfile for D escape flux
    p: path in which the results are located
    titleext: optional string to append to the plot title
    poster: whether to make font big or not
    =#

    pth = p*"/"
    const times = h5read(rfH,"fluxes/times")  # read in the times
    const fluxH = h5read(rfH,"fluxes/fluxvals")  # read in the flux values for H
    const fluxD = h5read(rfD,"fluxes/fluxvals")  # read in the flux values

    # for each ppm, fill dict with 1D array of variation vs altitude (altitudes are the dict keys)
    flux80H = Dict{Float64,Array{Float64, 1}}()
    flux20H = Dict{Float64,Array{Float64, 1}}()
    flux40H = Dict{Float64,Array{Float64, 1}}()
    flux60H = Dict{Float64,Array{Float64, 1}}()
    flux80D = Dict{Float64,Array{Float64, 1}}()
    flux20D = Dict{Float64,Array{Float64, 1}}()
    flux40D = Dict{Float64,Array{Float64, 1}}()
    flux60D = Dict{Float64,Array{Float64, 1}}()

    # also, the concentrations of H2O and HDO at the ground. each is a time series.
    h2o80 = Dict{Float64,Array{Float64, 1}}()
    h2o20 = Dict{Float64,Array{Float64, 1}}()
    h2o40 = Dict{Float64,Array{Float64, 1}}()
    h2o60 = Dict{Float64,Array{Float64, 1}}()
    hdo80 = Dict{Float64,Array{Float64, 1}}()
    hdo20 = Dict{Float64,Array{Float64, 1}}()
    hdo40 = Dict{Float64,Array{Float64, 1}}()
    hdo60 = Dict{Float64,Array{Float64, 1}}()

    for col_i in range(1, 6)
        # in these pairs, each dictionary key is an col_i (e.g. 40) and each value is the
        # list of ppms. index order is sheet, column, row. sht 1 80ppm; sht 2 20ppm;
        # sht 3 40 ppm; sheet 4 60ppm

        # get the altitude for ease of use.
        alt = fluxH[1, col_i, 2]

        # ppm 80
        flux80H[alt] =  fluxH[1, col_i, :][3:end-1]
        flux80D[alt] =  fluxD[1, col_i, :][3:end-1]
        fn = pth*"ppm_80_alt_"*string(Int(alt))*".h5"
        h2o80[alt], hdo80[alt] = conc_over_time_water(fn)

        # ppm 20
        flux20H[alt] = fluxH[2, col_i, :][3:end-1]     # first argument should be 2 if doing all ppm
        flux20D[alt] = fluxD[2, col_i, :][3:end-1]     # first argument should be 2 if doing all ppm
        fn = pth*"ppm_20_alt_"*string(Int(alt))*".h5"
        h2o20[alt], hdo20[alt] = conc_over_time_water(fn)

        # ppm 40
        flux40H[alt] = fluxH[3, col_i, :][3:end-1]
        flux40D[alt] = fluxD[3, col_i, :][3:end-1]
        fn = pth*"ppm_40_alt_"*string(Int(alt))*".h5"
        h2o40[alt], hdo40[alt] = conc_over_time_water(fn)

        # ppm 60
        flux60H[alt] = fluxH[4, col_i, :][3:end-1]
        flux60D[alt] = fluxD[4, col_i, :][3:end-1]
        fn = pth*"ppm_60_alt_"*string(Int(alt))*".h5"
        h2o60[alt], hdo60[alt] = conc_over_time_water(fn)
    end

    # make some arrays to hold all the R values --------------------------------
    R20_all = Float64[]
    R40_all = Float64[]
    R60_all = Float64[]
    R80_all = Float64[]

    alt_all = Array{Float64}(6)

    # stuff for all plots ------------------------------------------------------
    keyz = sort(collect(keys(flux20H)))
    colors = ["cornflowerblue", "blue", "silver", "dimgray", "salmon", "red"]
    lbls = ["20 km", "40 km", "60 km", "80 km", "100 km", "120 km"]
    linestyles = ["-", "-", "-", "-"] #"--", "-.", ":"] # toggle for styles

    # PLOT ---------------------------------------------------------------------
    fig, ax = subplots(1, 2, figsize=(20,10))
    i = 1
    for k in keyz
        col = colors[i]
        lbl = lbls[i]

        R20 = 2*(flux20D[k]./flux20H[k])./(hdo20[k]./h2o20[k])
        R40 = 2*(flux40D[k]./flux40H[k])./(hdo40[k]./h2o40[k])
        R60 = 2*(flux60D[k]./flux60H[k])./(hdo60[k]./h2o60[k])
        R80 = 2*(flux80D[k]./flux80H[k])./(hdo80[k]./h2o80[k])

        # arrays of all the values for each ppm
        append!(R20_all, R20)
        append!(R40_all, R40)
        append!(R60_all, R60)
        append!(R80_all, R80)

        # averages for each altitude.
        alt_all[Int(k/20)] = mean(mean([R20, R40, R60, R80])) # mean first across ppm, then over time.

        # the short term plot
        ax[1, 1][:semilogx](times, R20, color=col, label=lbl, linestyle=linestyles[1])
        ax[1, 1][:semilogx](times, R40, color=col, linestyle=linestyles[2])
        ax[1, 1][:semilogx](times, R60, color=col, linestyle=linestyles[3])
        ax[1, 1][:semilogx](times, R80, color=col, linestyle=linestyles[4])

        # the long term plot
        ax[2, 1][:semilogx](times, R20, color=col, label=lbl, linestyle=linestyles[1])
        ax[2, 1][:semilogx](times, R40, color=col, linestyle=linestyles[2])
        ax[2, 1][:semilogx](times, R60, color=col, linestyle=linestyles[3])
        ax[2, 1][:semilogx](times, R80, color=col, linestyle=linestyles[4])
        i += 1
    end

    # control font sizes for poster making
    if poster==true
        fs = Dict("ticks"=>24, "labels"=>28, "legend"=>18, "title"=>30, "stitle"=>34)
    elseif poster==false
        fs = Dict("ticks"=>16, "labels"=>20, "legend"=>14, "title"=>22, "stitle"=>24)
    end

    # plot label stuff
    ax[1, 1][:tick_params](axis="both", which="both", labelsize=fs["ticks"])
    ax[1, 1][:set_xlabel]("Time (seconds)", fontsize=fs["labels"])
    ax[1, 1][:set_ylabel]("R value", fontsize=fs["labels"])
    ax[1, 1][:set_xlim](10^1,10^7)
    ax[1, 1][:xaxis][:grid](which="major", color="gainsboro")
    ax[1, 1][:yaxis][:grid](which="minor", color="gainsboro")
    ax[1, 1][:set_title]("1 year", fontsize=fs["title"])
    ax[2, 1][:tick_params](axis="both", which="both", labelsize=fs["ticks"])
    ax[2, 1][:set_xlabel]("Time (seconds)", fontsize=fs["labels"])
    ax[2, 1][:set_ylabel]("R value", fontsize=fs["labels"])
    ax[2, 1][:set_xlim](10^1,10^15)
    ax[2, 1][:xaxis][:grid](which="major", color="gainsboro")
    ax[2, 1][:yaxis][:grid](which="minor", color="gainsboro")
    ax[2, 1][:set_title]("10 My", fontsize=fs["title"])
    legend(fontsize=fs["legend"], loc="upper left")

    tight_layout(rect=[0, 0.03, 1, 0.95]) # use if you want to do a suptitle: rect=[0, 0.03, 1, 0.95]

    # get the temperatures
    junk = split(titleext, "_")
    if junk[1]=="temp"
        suptitle("D/H Fractionation, "*L"T_{surf}="*"$(junk[2]), "*L"T_{tropo}="*"$(junk[3]), "*L"T_{exo}="*"$(junk[4])", fontsize=fs["stitle"])
    elseif junk[1]=="water"
        suptitle("D/H Fractionation, "*L"Water fraction="*"$(junk[2])", fontsize=fs["stitle"])
    else
        println("PROBLEM")
    end

    fnbase = "R"
    if poster==true
        savefig(pth*fnbase*"_poster.png", bbox_inches="tight")
    elseif poster==false
        savefig(pth*fnbase*".png", bbox_inches="tight")
    end

    # R averages --------------------------------------------------------------
    open(pth*"R_averages.txt", "w") do f
        write(f, "Avg R for 20 ppm: $(@sprintf("%.2e", mean(R20_all)))\n")
        write(f, "Avg R for 40 ppm: $(@sprintf("%.2e", mean(R40_all)))\n")
        write(f, "Avg R for 60 ppm: $(@sprintf("%.2e", mean(R60_all)))\n")
        write(f, "Avg R for 80 ppm: $(@sprintf("%.2e", mean(R80_all)))\n")

        for i in range(1,6)
            write(f, "Avg for alt $(i*20) km: $(@sprintf("%.2e", alt_all[i]))\n")
        end
    end
end

function calcFractionation_single(rfH, rfD, p, caseno="", titleext="", poster=false)
    pth = p*"Case "*caseno*"/"
    const times = h5read(rfH,"fluxes/times")  # read in the times
    const fluxH_data = h5read(rfH,"fluxes/fluxvals")  # read in the flux values for H
    const fluxD_data = h5read(rfD,"fluxes/fluxvals")  # read in the flux values

    # for each ppm, fill dict with 1D array of variation vs altitude (altitudes are the dict keys)
    fluxH = Dict{Float64,Array{Float64, 1}}()
    fluxD = Dict{Float64,Array{Float64, 1}}()

    # also, the concentrations of H2O and HDO at the ground. each is a time series.
    h2opop = Dict{Float64,Array{Float64, 1}}()
    hdopop = Dict{Float64,Array{Float64, 1}}()

    # get the altitude for ease of use.
    alt = fluxH_data[1, 1, 2]
    ppm = fluxH_data[1, 1, 1]

    # get the fluxes and the water concentrations at the ground
    fluxH =  fluxH_data[1, 1, :][3:end-1]
    fluxD =  fluxD_data[1, 1, :][3:end-1]
    fn = pth*"ppm_"*string(Int(ppm))*"_alt_"*string(Int(alt))*".h5"
    h2opop, hdopop = conc_over_time_water(fn)

    # stuff or all plots -------------------------------------------------------
    # calculate fractionation factor
    ff = 2*(fluxD./fluxH)./(hdopop./h2opop)
    ff_meanval = mean(ff)
    println(ff_meanval)

    # control font sizes for poster making
    if poster==true
        fs = Dict("ticks"=>24, "labels"=>28, "legend"=>18, "title"=>30)
    elseif poster==false
        fs = Dict("ticks"=>16, "labels"=>20, "legend"=>14, "title"=>22)
    end

    # short term plot (1 year) ------------------------------------------------
    fig, ax = subplots(figsize=(12,8))
    ax[:semilogx](times, ff, color="blue")
    xlabel("Time (seconds)", fontsize=fs["labels"])
    ylabel("R value", fontsize=fs["labels"])
    xlim(10^0,10^7)
    #ylim(ff_meanval - 0.1*ff_meanval, ff_meanval + 0.1*ff_meanval, )
    xticks(fontsize = fs["ticks"])
    yticks(fontsize = fs["ticks"])
    ax[:xaxis][:grid](which="major", color="gainsboro")
    ax[:yaxis][:grid](which="minor", color="gainsboro")
    title("Fractionation factor "*titleext, fontsize=fs["title"])

    fnbase = "R"
    if poster==true
        savefig(pth*fnbase*"_"*"Case"*caseno*"_poster.png")
    elseif poster==false
        savefig(pth*fnbase*"_"*"Case"*caseno*".png")
    end

    # long term plot -----------------------------------------------------------
    fig, ax = subplots(figsize=(12,8))
    ax[:semilogx](times, ff, color="blue")
    xlabel("Time (seconds)", fontsize=fs["labels"])
    ylabel("Fractionation factor", fontsize=fs["labels"])
    xlim(10^1, 10^17)
    # ylim(ff_meanval - 0.1*ff_meanval, ff_meanval + 0.1*ff_meanval, )
    xticks([10^3, 10^6, 10^9, 10^12, 10^15], fontsize=fs["ticks"])
    yticks(fontsize=fs["ticks"])
    ax[:xaxis][:grid](which="major", color="gainsboro")
    ax[:yaxis][:grid](which="minor", color="gainsboro")
    title("Fractionation Factor "*titleext, fontsize=fs["title"])

    fnbase = "R_long"
    if poster==true
        savefig(pth*fnbase*"_"*"Case"*caseno*"_poster.png")
    elseif poster==false
        savefig(pth*fnbase*"_"*"Case"*caseno*".png")
    end
end

function calcPopRatio(p)
    #=
    Calculates and plots the ratio (HD/H2)/(HDO/H2O)
    =#

    # *** TODO: FIX THIS FUNCTION IT IS NOT WORKING BECAUSE MESSY PASTING
    # set up the filenames
    waterppmvec = [20 40 60 80]
    wateraltvec = [20 40 60 80 100 120]
    parmsvec = [[a,b] for a in waterppmvec, b in wateraltvec]
    parmsvec = reshape(parmsvec,length(parmsvec))
    files = [string(p*"ppm_",a[1],"_alt_",a[2],".h5") for a in parmsvec]
    const times = h5read(files[1],"n_current/timelist")

    println("Ending the setup: ", Dates.Time(now()))

    println("Starting the zip: ", Dates.Time(now()))
    popratios = Dict{String, Array{Float64, 2}}()

    # loop through files collecting concentrations over time
    for (pair, f) in zip(parmsvec, files)
        # get the concentrations of H2O, HDO, H2, HD at ground level
        concarray = Array{Float64}(4,1299)
        concarray[1,:] = conc_over_time_water(f)[1] # H2O
        concarray[2,:] = conc_over_time_water(f)[2] # HDO
        concarray[3,:] = conc_over_time_molecH(f)[1] # H2
        concarray[4,:] = conc_over_time_molecH(f)[2] # HD
        popratios[string(pair[1])*"_"*string(pair[2])] = concarray
    end


    # stuff for all plots ------------------------------------------------------
    colors = ["cornflowerblue", "blue", "silver", "dimgray", "salmon", "red"]
    lbls = ["20 km", "40 km", "60 km", "80 km", "100 km", "120 km"]
    linestyles = ["-", "--", "-.", ":"]# "-", "-", "-"] #toggle for styles

    means_by_alt = Array{Float64}(6, 1)

    # plotting labels and things
    if poster==true
        fs = Dict("ticks"=>24, "labels"=>28, "legend"=>18, "title"=>30)
    elseif poster==false
        fs = Dict("ticks"=>16, "labels"=>20, "legend"=>14, "title"=>22)
    end

    fig, ax = subplots(figsize=(12,8))
    title(latexstring("(HD/H_2) / (HDO/H_2O)"), fontsize=fs["title"])
    xticks(fontsize = fs["ticks"])
    yticks(fontsize = fs["ticks"])
    ax[:xaxis][:grid](which="major", color="gainsboro")
    ax[:yaxis][:grid](which="both", color="gainsboro")
    xlabel("Time (seconds)", fontsize=fs["labels"])
    ylabel("Ratio", fontsize=fs["labels"])


    # short term plot (1 year)
    i = 1
    for a in wateraltvec  # loop over altitudes
        thisalt_ratios = Array{Float64}(4,1299)
        k = 1
        for ppm in waterppmvec  # loop over ppm

            keystr = string(ppm)*"_"*string(a)
            conc = popratios[keystr]  # order is H2O, HDO, H2, HD

            col = colors[i]
            lbl = lbls[i]

            # calculate the population ratio for this altitude and store it for averaging
            thisratio = (conc[4, :]./conc[3, :]) ./ (conc[2, :]./conc[1, :])
            thisalt_ratios[k, :] = thisratio

            # plot
            ax[:semilogx](times, thisratio, color=col, label=lbl, linestyle=linestyles[k])

            k += 1
        end
        # averages for each altitude.
        means_by_alt[i] = mean(thisalt_ratios) # mean across ppm and time.

        i += 1 # go to the next color and label (based on altitudes)

    end

    println(means_by_alt)

    show()
end

lead = "/data/GoogleDrive/"#"/home/emc/GoogleDrive/"#

# for doing Yung experiments
# P = lead*"Phys/LASP/chaffincode-working/Yung-With-Old-Water/"#
# println("Input case number (just the number): ")
# caseno = readline()#"Case 1"#
# fH = P*"Case "*caseno*"/H_esc_flux_history.h5"
# fD = P*"Case "*caseno*"/D_esc_flux_history.h5"
# calcFractionation_single(fH, fD, P, caseno)

# for our experiments
expfolder = isdefined(:ARGS) ? ARGS[1] : println("Please use command line args")
P = "/data/VaryTW_Ana/"*expfolder
fH = P*"/H_esc_flux_history.h5"
fD = P*"/D_esc_flux_history.h5"
calcFractionation(fH, fD, P, expfolder, true)

#calcPopRatio(P)
