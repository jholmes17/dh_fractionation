using PyPlot
using HDF5
using LaTeXStrings

function plot_Hflux(pth, readfile, dincluded, expfolder, poster=false)
    #=
    Creates a plot of H flux out of the atmosphere over time
    pth: parent folder of readfile. used to save figures properly
    readfile: file of H fluxes out of the atmosphere
    dincluded: whether it includes D fluxes. values: "H" or "H+D"
    poster: whether to use text sizes appropriate for a poster plot
    =#

    const times = h5read(readfile,"fluxes/times")  # read in the times
    const flux = h5read(readfile,"fluxes/fluxvals")  # read in the flux values

    # for each center of the Gaussian ppm distribution, fill dict with 1D array of variation
    # vs altitude (altitudes are the dict keys)
    flux80 = Dict{Float64,Array{Float64, 1}}()
    flux20 = Dict{Float64,Array{Float64, 1}}()
    flux40 = Dict{Float64,Array{Float64, 1}}()
    flux60 = Dict{Float64,Array{Float64, 1}}()

    for altitude in range(1, 6)
        # in these pairs, each dictionary key is an altitude (e.g. 40) and each value is the
        # list of ppms. index order is sheet, column, row. sht 1 80ppm; sht 2 20ppm;
        # sht 3 40 ppm; sheet 4 60ppm

        a = flux[4, altitude, :][2]
        flux80[a] = flux[1, altitude, :][3:end-1]
        flux20[a] = flux[2, altitude, :][3:end-1]  # !! CHANGE THE 1 BACK TO "2" WHEN DOING ALL PROFILES
        flux40[a] = flux[3, altitude, :][3:end-1]
        flux60[a] = flux[4, altitude, :][3:end-1]
    end

    # stuff for all plots
    keyz = sort(collect(keys(flux20)))
    colors = ["cornflowerblue", "blue", "silver", "dimgray", "salmon", "red"]
    lbls = ["20 km", "40 km", "60 km", "80 km", "100 km", "120 km"]
    linestyles = ["-", "-", "-", "-"]

    # short term plot (1 year) -------------------------------------------------
    fig, ax = subplots(figsize=(12,8))
    i = 1
    miny = 1e8
    maxy = 0 # find the maxy to set plot limits
    for k in keyz
        col = colors[i]
        lbl = lbls[i]# * " " *
        ax[:loglog](times, flux20[k], color=col, label=lbl,
                    linestyle=linestyles[1])
        ax[:loglog](times, flux40[k], color=col, linestyle=linestyles[2])
        ax[:loglog](times, flux60[k], color=col, linestyle=linestyles[3])
        ax[:loglog](times, flux80[k], color=col, linestyle=linestyles[4])

        # find the highest and lowest values for setting plot limits
        curmax = maximum([maximum(flux20[k]), maximum(flux40[k]), maximum(flux60[k]), maximum(flux80[k])])
        curmin = minimum([minimum(flux20[k]), minimum(flux40[k]), minimum(flux60[k]), minimum(flux80[k])])
        if curmin < miny
            miny = 10^floor(log10(curmin))
        end
        if curmax > maxy
            maxy = 10^ceil(log10(curmax))
        end
        i += 1
    end
    if miny == maxy
        maxy *= 10
    end

    # control font sizes for poster making
    if poster==true
        fs = Dict("ticks"=>24, "labels"=>28, "legend"=>18, "title"=>30, "stitle"=>34)
    elseif poster==false
        fs = Dict("ticks"=>16, "labels"=>20, "legend"=>14, "title"=>22, "stitle"=>24)
    end

    # get the temperatures
    junk = split(expfolder, "_")
    if junk[1]=="temp"
        titleext = L"T_{surf}="*"$(junk[2]), "*L"T_{tropo}="*"$(junk[3]), "*L"T_{exo}="*"$(junk[4])"
    elseif junk[1]=="water"
        titleext = L"Water fraction="*"$(junk[2])"
    else
        println("PROBLEM")
    end

    # labels and stuff
    tight_layout(rect=[0, 0.03, 1, 0.95]) # use if you want to do a suptitle: rect=[0, 0.03, 1, 0.95]
    xlabel("Time (seconds)", fontsize=fs["labels"])
    ylabel(dincluded*" escape flux ("*latexstring("cm^{-2} s^{-1}")*")",
           fontsize=fs["labels"])
    xlim(10^3,10^7)
    ylim(miny, maxy)
    xticks(fontsize = fs["ticks"])
    yticks(fontsize = fs["ticks"])
    ax[:xaxis][:grid](which="major", color="gainsboro")
    ax[:yaxis][:grid](which="minor", color="gainsboro")
    legend(fontsize=fs["legend"], loc="upper left")
    suptitle("Flux response to elevated "*latexstring("H_2O, HDO"), fontsize=fs["stitle"], y=1.05)
    title(titleext, fontsize=fs["title"])
    if poster==true
        savefig(pth*"atmoresponse" * "_" * dincluded * "_poster.png", bbox_inches="tight")
    elseif poster==false
        savefig(pth*"atmoresponse" * "_" * dincluded * ".png", bbox_inches="tight")
    end

    # long term plot -----------------------------------------------------------
    fig, ax = subplots(figsize=(12,8))
    i = 1
    for k in keyz
        col = colors[i]
        lbl = lbls[i]
        ax[:loglog](times, flux20[k], color=col, label=lbl,
                    linestyle=linestyles[1])
        ax[:loglog](times, flux40[k], color=col, linestyle=linestyles[2])
        ax[:loglog](times, flux60[k], color=col, linestyle=linestyles[3])
        ax[:loglog](times, flux80[k], color=col, linestyle=linestyles[4])

        # find the highest and lowest values for setting plot limits
        curmax = maximum([maximum(flux20[k]), maximum(flux40[k]), maximum(flux60[k]), maximum(flux80[k])])
        curmin = minimum([minimum(flux20[k]), minimum(flux40[k]), minimum(flux60[k]), minimum(flux80[k])])
        if curmin < miny
            miny = 10^floor(log10(curmin))
        end
        if curmax > maxy
            maxy = 10^ceil(log10(curmax))
        end
        i += 1
    end
    if miny == maxy
        maxy *= 10
    end

    # get the temperatures
    junk = split(expfolder, "_")
    if junk[1]=="temp"
        titleext = L"T_{surf}="*"$(junk[2]), "*L"T_{tropo}="*"$(junk[3]), "*L"T_{exo}="*"$(junk[4])"
    elseif junk[1]=="water"
        titleext = L"Water fraction="*"$(junk[2])"
    else
        println("PROBLEM")
    end

    # labels and stuff
    tight_layout(rect=[0, 0.03, 1, 0.95]) # use if you want to do a suptitle: rect=[0, 0.03, 1, 0.95]
    xlabel("Time (seconds)", fontsize=fs["labels"])
    ylabel(dincluded*" escape flux ("*latexstring("cm^{-2} s^{-1}")*")",
           fontsize=fs["labels"])
    xlim(10^1, 10^17)
    ylim(miny, maxy)
    xticks([10^3, 10^6, 10^9, 10^12, 10^15], fontsize=fs["ticks"])
    yticks(fontsize=fs["ticks"])
    ax[:xaxis][:grid](which="major", color="gainsboro")
    ax[:yaxis][:grid](which="minor", color="gainsboro")
    legend(fontsize=fs["legend"], loc="upper left")
    suptitle("Flux response to elevated "*latexstring("H_2O, HDO"), fontsize=fs["stitle"], y=1.05)
    title(titleext, fontsize=fs["title"])
    if poster==true
        savefig(pth*"atmoresponse" * "_" * dincluded * "_long_poster.png", bbox_inches="tight")
    elseif poster==false
        savefig(pth*"atmoresponse" * "_" * dincluded * "_long.png", bbox_inches="tight")
    end
end

lead = "/data/VaryTW_Ana/"#"/home/emc/GoogleDrive/"#
expfolder = isdefined(:ARGS) ? ARGS[1] : println("Please use command line args")
fbase = lead*expfolder*"/"
fH = fbase*"H_esc_flux_history.h5"
plot_Hflux(fbase, fH, "H", expfolder, true)
fHD = fbase*"H_and_D_esc_flux_history.h5"
plot_Hflux(fbase, fHD, "H+D", expfolder, true)
fD = fbase*"D_esc_flux_history.h5"
plot_Hflux(fbase, fD, "D", expfolder, true)


# for replotting Mike's work
# mikefile = "/home/emc/GoogleDrive/Phys/LASP/Mikes results/outputfiles_version06_original/H_esc_flux_history.h5"
# plot_Hflux(mikefile, "H", false)
