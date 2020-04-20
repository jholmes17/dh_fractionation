################################################################################
# plot_water_loss.jl
# TYPE: Analysis
# WHICH: Equilibrium
# DESCRIPTION: Estimates water lost from Mars for a given array of fractionation 
# factors and plots it, along with the Rayleigh distillation equation.
#
# Updated 26 December 2019
# Eryn Cangi
# Currently tested for Julia: 0.7
################################################################################

using Analysis
using PyPlot
using PyCall

include("PARAMETERS.jl")

function h2olost(cur_h2o, DHnow=5.5, DH0=1, f=0.001)
    return cur_h2o * ((DHnow/DH0)^(1/(1-f)) - 1)
end

function required_escape_rate(W_lost)
    #=
    Calculates the required escape rate of H atoms, assuming a constant escape
    rate and a linear change with time. NOTE: this assumes all water lost 
    (W_lost) is in the form H2O.
    =#
    H_lost = GEL_to_molecule(W_lost, "H")
    flux = H_lost / 1.419e17 # units #/cm^2s 
end

function plot_water_loss(f_therm, f_both, f_therm_range=nothing, f_therm_both=nothing)
    #=
    Plot water loss as a result of fractionation factor results and also a line
    according to the D/H ratio used. 

    f_therm: mean fractionation factor for only thermal escape
    f_both: mean fractionation factor for both thermal + nonthermal
    f_therm_range: min and max fractionation factor for thermal escape
    f_both_range: min and max fractionation factor for all types of escape
    =#
    # set up the domain 
    cur_h2o = range(20, stop=30, length=11)

    # make plot look nice
    fig = figure(figsize=(6,4))
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = "Louis George Caf?"
    rcParams["font.monospace"] = "FreeMono"
    rcParams["font.size"] = 26
    rcParams["axes.labelsize"] = 28
    rcParams["xtick.labelsize"] = 26
    rcParams["ytick.labelsize"] = 26
    ax = gca()
    plot_bg(ax)
    xlabel("Current exchangeable water (m GEL)")
    ylabel("Water lost to space (m GEL)")
    title("Water lost since 4.5 Ga")

    # plot 
    blues = get_grad_colors(2, "Blues", strt=0.5)

    # calculate and plot the loss for mean f
    loss_t = h2olost(cur_h2o, 5.5, 1.275, f_therm)
    loss_b = h2olost(cur_h2o, 5.5, 1.275, f_both)

    ax.plot(cur_h2o, loss_t, zorder=2, color=blues[1, :], marker="D", ms=12)
    ax.plot(cur_h2o, loss_b, zorder=2, color=blues[2, :], marker="D", ms=12)

    if f_therm_range != nothing
        therm_min_loss = h2olost(cur_h2o, 5.5, 1.275, f_therm_range[1])
        therm_max_loss = h2olost(cur_h2o, 5.5, 1.275, f_therm_range[2])
        ax.fill_between(cur_h2o, therm_min_loss, therm_max_loss, zorder=2,
                        alpha=0.4, color=blues[1, :])
    end

    if f_both_range != nothing
        both_min_loss = h2olost(cur_h2o, 5.5, 1.275, f_both_range[1])
        both_max_loss = h2olost(cur_h2o, 5.5, 1.275, f_both_range[2])
        ax.fill_between(cur_h2o, both_min_loss, both_max_loss, zorder=2,
                        alpha=0.4, color=blues[2, :])
    end
    println("Thermal escape water loss (based on f for mean temperature profile):")
    println(Array(loss_t))
    println("All escape water loss (based on f for mean temperature profile): ")
    println(Array(loss_b))

    println()

    println("Thermal escape minimum water loss (f minima):")
    println(Array(therm_min_loss))
    println("Thermal escape maximum water loss (f maxima):")
    println(Array(therm_max_loss))

    println()

    println("All escape minimum water loss (f minima):")
    println(Array(both_min_loss))
    println("All escape maximum water loss (f maxima):")
    println(Array(both_max_loss))

    text(21, 83, "thermal+non-thermal escape", rotation=24, color=blues[2, :])
    text(22, 68, "thermal escape only", rotation=23, color=blues[1, :])

    savefig("../Results/ALL STUDY PLOTS/h2oloss.png", bbox_inches="tight")
end

function plot_simple_water_loss(cur_h2o=25)
    #=
    Plot just the Rayleigh distillation equation, no fractionation factor results.
    For demonstration purposes.

    cur_h2o: assumed water inventory in m GEL.
    =#
    farray = range(0.000005, stop=1, length=100)

    l1 = [h2olost(cur_h2o, 5.5, 1.275, f) for f in farray]

    # make plot
    fig = figure(figsize=(10, 7))

    # style.use("default")
    rc("text", usetex=false)
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = "Louis George Caf?"
    rcParams["font.monospace"] = "FreeMono"
    rcParams["font.size"] = 26
    rcParams["axes.labelsize"] = 28
    rcParams["xtick.labelsize"] = 26
    rcParams["ytick.labelsize"] = 26

    ax = gca()
    plot_bg(ax)

    xlabel(r"Fractionation factor $f$")
    ylabel("Watr lost (m GEL) in 4.5 Ga")

    # plot things
    ax.semilogx(farray, l1, zorder=2, color="cornflowerblue")
    # legend(fontsize=24)
    text(0.000005, 50, "Water loss constant \n"+L" with $f$ when $f\ll0.05$")
    text(0.0006, 200, "Water loss depends \nstrongly on " + L"$f$ for $f\gtrapprox 0.05$")
    ylim(0, 240)
    savefig("../Results/ALL STUDY PLOTS/h2oloss_simple.png", bbox_inches="tight")
end

function main()
    println("ALERT: Make sure the f values are correct."*
            " They must be filled in manually, using output from plot_f_results.jl")
    f_mean_thermal = 0.0018797755997056428
    f_therm_range = [3.23715e-5, 0.0170131]
    f_mean_both = 0.05971712335227322
    f_both_range = [0.032912, 0.102279]

    plot_water_loss(f_mean_thermal, f_mean_both, f_therm_range, f_both_range)
end

