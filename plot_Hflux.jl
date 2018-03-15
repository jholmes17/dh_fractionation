using HDF5
using PyPlot
using LaTeXStrings

function plot_Hflux(readfile, dincluded, poster=false)
    #=
    Creates a plot of H flux out of the atmosphere over time
    readfile: file of H fluxes out of the atmosphere
    packettype: standard or high water
    dincluded: whether it includes D fluxes. values: "H" or "H+D"
    =#

    const times = h5read(readfile,"fluxes/times")  # read in the times
    const flux = h5read(readfile,"fluxes/fluxvals")  # read in the flux values

    # for each center of the Gaussian ppm distribution, fill dict with 1D array of variation
    # vs altitude (altitudes are the dict keys)
    ppm80 = Dict{Float64,Array{Float64, 1}}()
    ppm20 = Dict{Float64,Array{Float64, 1}}()
    ppm40 = Dict{Float64,Array{Float64, 1}}()
    ppm60 = Dict{Float64,Array{Float64, 1}}()

    for altitude in range(1, 6)
        # in these pairs, each dictionary key is an altitude (e.g. 40) and each value is the
        # list of ppms. index order is sheet, column, row. sht 1 80ppm; sht 2 20ppm;
        # sht 3 40 ppm; sheet 4 60ppm
        entries80 = flux[1, altitude, :]
        ppm80[entries80[2]] =  entries80[3:end-1]

        entries20 = flux[2, altitude, :]
        ppm20[entries20[2]] = entries20[3:end-1]

        entries40 = flux[3, altitude, :]
        ppm40[entries40[2]] = entries40[3:end-1]

        entries60 = flux[4, altitude, :]
        ppm60[entries60[2]] = entries60[3:end-1]
    end

    # stuff for all plots
    keyz = sort(collect(keys(ppm20)))
    colors = ["cornflowerblue", "blue", "silver", "dimgray", "salmon", "red"]
    lbls = ["20 km", "40 km", "60 km", "80 km", "100 km", "120 km"]
    linestyles = ["-", "-", "-", "-"]

    # short term plot (1 year)
    fig, ax = subplots(figsize=(12,8))
    i = 1
    for key in keyz
        col = colors[i]
        lbl = lbls[i]# * " " *
        ax[:loglog](times, ppm20[key], color=col, label=lbl,
                    linestyle=linestyles[1])
        ax[:loglog](times, ppm40[key], color=col, linestyle=linestyles[2])
        ax[:loglog](times, ppm60[key], color=col, linestyle=linestyles[3])
        ax[:loglog](times, ppm80[key], color=col, linestyle=linestyles[4])
        i += 1
    end

    if poster==true
        fs = Dict("ticks"=>24, "labels"=>28, "legend"=>18, "title"=>30)
    elseif poster==false
        fs = Dict("ticks"=>16, "labels"=>20, "legend"=>14, "title"=>22)
    end

    xlabel("Time (seconds)", fontsize=fs["labels"])
    ylabel(dincluded*" escape flux ("*latexstring("cm^{-2} s^{-1}")*")",
           fontsize=fs["labels"])
    xlim(10^3,10^7)
    if dincluded=="D"
        ylim(10^2,10^5)
    else
        ylim(10^8, 10^10)
    end

    xticks(fontsize = fs["ticks"])
    yticks(fontsize = fs["ticks"])

    ax[:xaxis][:grid](which="major", color="gainsboro")
    ax[:yaxis][:grid](which="minor", color="gainsboro")
    legend(fontsize=fs["legend"], loc="upper left")
    title("Atmospheric response to elevated "*latexstring("H_2O, HDO"), fontsize=fs["title"])
    if poster==true
        savefig("atmoresponse" * "_" * dincluded * "_poster.png")
    elseif poster==false
        savefig("atmoresponse" * "_" * dincluded * ".png")
    end

    #long term plot
    fig, ax = subplots(figsize=(12,8))
    keyz = sort(collect(keys(ppm20)))
    i = 1
    for key in keyz
        col = colors[i]
        lbl = lbls[i]
        ax[:loglog](times, ppm20[key], color=col, label=lbl,
                    linestyle=linestyles[1])
        ax[:loglog](times, ppm40[key], color=col, linestyle=linestyles[2])
        ax[:loglog](times, ppm60[key], color=col, linestyle=linestyles[3])
        ax[:loglog](times, ppm80[key], color=col, linestyle=linestyles[4])
        i += 1
    end
    xlabel("Time (seconds)", fontsize=fs["labels"])
    ylabel(dincluded*" escape flux ("*latexstring("cm^{-2} s^{-1}")*")",
           fontsize=fs["labels"])
    xlim(10^1, 10^17)
    if dincluded=="D"
        ylim(10^2,10^5)
    else
        ylim(10^8, 10^10)
    end
    xticks([10^3, 10^6, 10^9, 10^12, 10^15], fontsize=fs["ticks"])
    yticks(fontsize=fs["ticks"])

    ax[:xaxis][:grid](which="major", color="gainsboro")
    ax[:yaxis][:grid](which="minor", color="gainsboro")
    legend(fontsize=fs["legend"], loc="upper left")
    title("Atmospheric response to elevated "*latexstring("H_2O, HDO"), fontsize=fs["title"])
    if poster==true
        savefig("atmoresponse" * "_" * dincluded * "_long_poster.png")
    elseif poster==false
        savefig("atmoresponse" * "_" * dincluded * "_long.png")
    end
end

#Results-Standard Water/
f1 = "/home/emc/Google Drive/Phys/LASP/Mars/chaffincode-working/H_esc_flux_history.h5"
plot_Hflux(f1, "H", true)
f2 = "/home/emc/Google Drive/Phys/LASP/Mars/chaffincode-working/H_and_D_esc_flux_history.h5"
plot_Hflux(f2, "H+D", true)
f3 = "/home/emc/Google Drive/Phys/LASP/Mars/chaffincode-working/D_esc_flux_history.h5"
plot_Hflux(f3, "D", true)

# f4 = "/home/emc/Google Drive/Phys/LASP/Mars/chaffincode v0.4.7/H_esc_flux_history-Mikesresult.h5"
# plot_Hflux(f4, "H", true)
