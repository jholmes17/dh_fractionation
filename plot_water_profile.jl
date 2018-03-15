using HDF5
using PyPlot

function plot_water_profile(readfile, packettype)
    #=
    Creates a plot of water concentration in PPM vs. altitude.
    =#

    const alt = h5read(readfile,"waterprofs/alt")  # read in the altitudes
    const ppm = h5read(readfile,"waterprofs/ppm")  # read in the ppm

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
        entries80 = ppm[1, altitude, :]
        ppm80[entries80[2]] =  entries80[3:end]

        entries20 = ppm[2, altitude, :]
        ppm20[entries20[2]] = entries20[3:end]

        entries40 = ppm[3, altitude, :]
        ppm40[entries40[2]] = entries40[3:end]

        entries60 = ppm[4, altitude, :]
        ppm60[entries60[2]] = entries60[3:end]
    end

    # plot the stuff
    y = alt./1E5
    fig, ax = subplots(figsize=(6,8))
    keyz = sort(collect(keys(ppm20)))
    colors = ["cornflowerblue", "blue", "silver", "dimgray", "salmon", "red"]
    i = 1
    for key in keyz
        color = colors[i]
        plot(ppm20[key], y, color)
        plot(ppm40[key], y, color)
        plot(ppm60[key], y, color)
        plot(ppm80[key], y, color)
        i += 1
    end

    # place nice labels of ppm on plot
    for ppm in collect([20,40,60,80])
        plot([ppm, ppm], [145, 155], "k-", lw=0.5)
        text(ppm-5, 160, "$(ppm)", fontsize=14)
    end
    text(7, 160, "+", fontsize=14)
    text(90,160,"ppm",fontsize=14)

    # place nice labels of altitude on plot
    for a in collect([20,40,60,80,100,120])
        if a==20
            x = 135
            y = 25
        else
            x = 100
            y = a
        end
        text(x, y, "$(a) km", color=colors[round(Int,a/20)], fontsize=14)
    end

    # manage axes/general plot features
    ylabel("Altitude (km)", fontsize=14)
    xlabel("Water concentration, ppmv", fontsize=14)
    xlim(0,200)
    ylim(0,200)
    xticks(fontsize=12)
    yticks(fontsize=12)
    title("Water profile assumed", fontsize=16)

    Mx = matplotlib[:ticker][:MultipleLocator](100) # Define interval of major ticks
    f = matplotlib[:ticker][:FormatStrFormatter]("%d") # Define format of tick labels
    ax[:xaxis][:set_major_locator](Mx) # Set interval of major ticks
    ax[:xaxis][:set_major_formatter](f) # Set format of tick labels
    savefig("waterprofile.png")
end

f1 = "/home/emc/Google Drive/Phys/LASP/Mars/chaffincode-working/Results - Standard Water/H_esc_flux_history.h5"
plot_water_profile(f1)
