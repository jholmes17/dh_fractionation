#=
Created 25 July 2018 to make plots for DPS 2018
Eryn Cangi
=#

using PyPlot
using HDF5
using LaTeXStrings

function plot_flux_and_ff(rootpath)
    #=
    From a selection of experiment folders, study the converged cases and
    report the fractionation factors and fluxes of H and D out of the atmosphere.

    rootpath: path containing many experiment folders
    =#

    # UTILITY FUNCTIONS AND CONSTANTS ==========================================
    # fundamental constants
    const boltzmannK = 1.38e-23;    # J/K
    const bigG = 6.67e-11;          # N m^2/kg^2
    const mH = 1.67e-27;            # kg
    const marsM = 0.1075*5.972e24;  # kg
    const radiusM = 3396e5;         # cm

    function get_pop(rf, alt, species)
        #=
        Get the population of a particular species at a particular altitude.
        rf: a file to read
        alt: altitude in km
        species: whatever species you want
        =#

        # identify the row number index for the species
        const slist = h5read(rf, "n_current/species") # get the list of species
        s_row_idx = findin(slist, [species])[1]

        # get the concentrations of the species at the specified level
        const concs = h5read(rf, "n_current/n_current_mat")

        # return statements depending on altitude
        if typeof(alt) == Int64
            return concs[alt/2 + 1, s_row_idx]
        elseif alt == "surface"
            return concs[1, s_row_idx]
        elseif alt == "exobase"
            return concs[end, s_row_idx]
        else
            println("Bad altitude selection! Please enter an altitude in km (multiples of 2):")
            alt = readline() / 2 + 1
        end
    end

    function effusion_velocity(Texo::Float64, m::Float64)
        #=
        Returns effusion velocity for a species.
        Texo: temperature of the exobase (upper boundary) in K
        m: mass of one molecule of species in amu
        =#
        lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+200e5))
        vth=sqrt(2*boltzmannK*Texo/(m*mH))
        v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)
        return v
    end

    function T_fromfile(z::Float64, tarray)
        #=
        Loads a prescribed temperature profile from a file. These files typically
        were collected from the MCD.
        z: altitude, in cm.
        tfile: file containing the temperature profile by altitude. must be in K,
               and the temp profile must be the second column of the file.
        =#
        idx = Int(z / 2e5) + 1
        return tarray[idx]
    end
    # END UTILITY FUNCTIONS AND CONSTANTS ======================================

    # get the list of convergence file paths using REGULAR EXPRESSIONS!! =======
    filelist = String[]
    explist = String[]
    temparrays = Array{Array{Float64}}(30)

    filepattern  = r"(converged_standardwater_D_).*(.h5)"
    exppattern = r"(?<=(D_)).*(?=(.h5))"

    i = 1
    for (root, dirs, files) in walkdir(rootpath)
        for file in files
            if ismatch(filepattern, file)
                expname = match(exppattern, file).match
                push!(filelist, rootpath*expname*"/"*file)
                push!(explist, expname)
                # load temp file cause we need it
                tfile = "/data/GoogleDrive/Phys/LASP/Resources/Water-Temp Profiles/"*expname*"-TEMP.dat"
                tfiledat = open(tfile, "r")
                tdata = readdlm(tfiledat, ' ', Float64)
                temparrays[i] = tdata[:, 2]
                i += 1
            end
        end
    end

    # calculate fluxes and fractionation factors ===============================
    fluxpairs = Array{Array{Float64}}(length(explist))  # a list of 2-el lists; each el will be [fluxD, fluxH]
    fracfacs = Array{Float64}(length(explist))          # just a list of numbers

    # set up the figure!
    fig, ax = subplots(1, 3, figsize=(15,5), sharex=true)#, sharey=true)
    ax1 = ax[1,1]
    ax2 = ax[2,1]
    ax3 = ax[3,1]
    ax1b = ax1[:twinx]()
    ax2b = ax2[:twinx]()
    ax3b = ax3[:twinx]()

    for i in 1:1:length(filelist)
        ex = explist[i]
        # load this experiment T profile and calculate effusion velocities
        Temp(z::Float64) = T_fromfile(z, temparrays[i])
        H_ev = effusion_velocity(Temp(196e5),1.0)
        H2_ev = effusion_velocity(Temp(196e5),2.0)
        D_ev = effusion_velocity(Temp(196e5),2.0)
        HD_ev = effusion_velocity(Temp(196e5),3.0)

        # get the populations at top of atmosphere (flux)
        fluxH = get_pop(filelist[i], "exobase", "H")*H_ev + 2*get_pop(filelist[i], "exobase", "H2")*H2_ev + get_pop(filelist[i], "exobase", "HD")*HD_ev
        fluxD = get_pop(filelist[i], "exobase", "HD")*HD_ev + get_pop(filelist[i], "exobase", "D")*D_ev
        groundH = get_pop(filelist[i], "surface", "H2O")
        groundD = get_pop(filelist[i], "surface", "HDO")
        fluxpairs[i] = [fluxD, fluxH]

        # calculate the fraction factor (flux ratio / pop ratio water)
        FF = (fluxD / fluxH) / (groundD / (2*groundH))
        fracfacs[i] = FF

        # get the minimum temp in the profile - will be at mid atmo
        T_mid = minimum([Temp(a) for a in 0:2e5:200e5])

        # set up the markers
        if contains(ex, "summer") && contains(ex, "morn")
            c = "maroon"
        elseif contains(ex, "summer") && contains(ex, "afternoon")
            c = "red"
        elseif contains(ex, "summer") && contains(ex, "night")
            c = "pink"
        elseif contains(ex, "winter") && contains(ex, "morn")
            c = "navy"
        elseif contains(ex, "winter") && contains(ex, "afternoon")
            c = "royalblue"
        elseif contains(ex, "winter") && contains(ex, "night")
            c = "lightskyblue"
        end
        if contains(ex, "eq")
            shp = "+"
        elseif contains(ex, "n-pole")
            shp = "*"
        elseif contains(ex, "n-midlat")
            shp = "^"
        elseif contains(ex, "s-midlat")
            shp = "v"
        elseif contains(ex, "s-pole")
            shp = "o"
        end

        # plot points on each plot: flux against surface, mid and exo temps
        ax1[:scatter](Temp(0.0), fluxD+fluxH, color="black", alpha=0.6, s=50, marker=shp, color=c)
        ax2[:scatter](T_mid, fluxD+fluxH, color="black", alpha=0.6, s=50, marker=shp, color=c)
        ax3[:scatter](Temp(196e5), fluxD+fluxH, color="black", alpha=0.6, s=50, marker=shp, color=c)

        ax1b[:scatter](Temp(0.0), FF, color="lightgray", alpha=0.6, s=50)
        ax2b[:scatter](T_mid, FF, color="lightgray", alpha=0.6, s=50)
        ax3b[:scatter](Temp(196e5), FF, color="lightgray", alpha=0.6, s=50)
    end
    # titles and labels
    ax1[:set_title]("Flux and Fractionation Factor vs. "*L"T_{surf}")
    ax2[:set_title]("Flux and Fractionation Factor vs. "*L"T_{tropopause}")
    ax3[:set_title]("Flux and Fractionation Factor vs. "*L"T_{exobase}")
    ax1[:set_ylabel]("H+D flux")
    ax3b[:set_ylabel]("Fractionation Factor "*L"(\phi_D/\phi_H) / ([HDO]/2[H_2O])", color="lightgray")

    # scaling and limits
    ax1[:set_yscale]("log")
    ax2[:set_yscale]("log")
    ax3[:set_yscale]("log")
    ax1b[:set_ylim](bottom=0)
    ax2b[:set_ylim](bottom=0)
    ax3b[:set_ylim](bottom=0)

    # ticks
    ax1b[:set_yticks]([])
    ax2b[:set_yticks]([])
    ax2[:set_yticks]([])
    ax3[:set_yticks]([])
    setp(ax3b[:get_yticklabels](),color="lightgray")
    my = matplotlib[:ticker][:MultipleLocator](0.005e8)
    ax1[:yaxis][:set_minor_locator](my)

    # aesthetics
    subplots_adjust(wspace=0.025, hspace=0)
    #rcParams.update({"font.size": 20})
    show()

end


plot_flux_and_ff("/data/VaryTW/")
