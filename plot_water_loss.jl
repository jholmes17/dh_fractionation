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

using Photochemistry
using PyPlot
using PyCall

function h2olost(cur_h2o, DHnow=5.5, DH0=1, f=0.001)
    return cur_h2o * ((DHnow/DH0)^(1/(1-f)) - 1)
end

function plot_water_loss(f_thermal, f_both; current_inventory=25, moredh=false, 
                         plot_past=false)
    #=
    Plot water loss as a result of fractionation factor results and also a line
    according to the D/H ratio used. 

    f_thermal: fractionation factor results for only thermal escape
    f_nonthermal: fractionation factor results for both thermal + nonthermal
    current_inventory: current water inventory in m GEL. If a vector, 
                        calculates lines for each.
    moredh: whether to plot a few other lines for different D/H values.
    plot_past: whether to plot past studies on the plot
    =#
    # set up the domain 
    farray = range(0.000005, stop=1, length=500)

    # make plot look nice
    fig = figure(figsize=(10, 7))
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
    xlabel(L"Fractionation factor $f$")
    ylabel("Water lost (m GEL)")
    title("Water lost since 4.5 Ga")

    # colors
    deemph = "#555555"  # color for de-emphasis

    # past studies
    if plot_past
        ax.scatter(0.32, 3.4, zorder=2, marker="v", s=100, color=deemph) # Yung
        ax.scatter(0.016, 65, zorder=2, marker="*", s=100, color=deemph) # Kras 2000 Model 1
        ax.scatter(0.135, 160, zorder=2,marker="^", s=100, color=deemph) # Kras 2000 Model 2
        text(0.1, 150, "K2000, Model 1", ha="right", color=deemph, fontsize=22)
        text(0.02, 60, "K2000, Model 2", color=deemph, fontsize=22)
        text(0.03, 0, "Yung 1988", color=deemph, fontsize=22)
    end

    if typeof(current_inventory)==Int64
        # case for a single water value 
        emph = "xkcd:cerulean blue"
        emphline = "xkcd:clear blue"
        loss = [h2olost(current_inventory, 5.5, 1.275, f) for f in farray]
        loss_thermal = [h2olost(current_inventory, 5.5, 1.275, f) for f in f_thermal]
        loss_both = [h2olost(current_inventory, 5.5, 1.275, f) for f in f_both]
        filename_tag = string(current_inventory)

        ax.semilogx(farray, loss, zorder=2, color=emphline)
        ax.scatter(f_thermal, loss_thermal, zorder=3, color=emph, marker="o", s=150,
                   label="Thermal escape only")
        ax.scatter(f_both, loss_both, zorder=3, color=emph, marker="D", s=150,
                   label="Thermal + non-thermal escape")
        ax.legend(fontsize=20)
    else
        # colors
        emph = ["xkcd:cerulean blue", "xkcd:vibrant blue", "xkcd:cobalt blue"]
        ls = [":", "--", "-"]

        # make the arrays to store lost water values
        loss = Array{Float64}(undef, length(current_inventory), 500)
        loss_thermal = Array{Float64}(undef, length(current_inventory), length(f_thermal))
        loss_both = Array{Float64}(undef, length(current_inventory), length(f_both))
        for i in range(1, stop=length(current_inventory))
            loss[i, :] = [h2olost(current_inventory[i], 5.5, 1.275, f) for f in farray]
            loss_thermal[i, :] = [h2olost(current_inventory[i], 5.5, 1.275, f) for f in f_thermal]
            loss_both[i, :] = [h2olost(current_inventory[i], 5.5, 1.275, f) for f in f_both]

            ax.semilogx(farray, loss[i, :], zorder=2, color=emph[i], linestyle=ls[i], 
                        label="Present-day water: $(current_inventory[i]) m GEL")
            ax.scatter(f_thermal, loss_thermal[i, :], zorder=3, color=emph[i], marker="o", s=150)
            ax.scatter(f_both, loss_both[i, :], zorder=3, color=emph[i], marker="D", s=150)
        end
        filename_tag = string(current_inventory[1]) * "-" * string(current_inventory[end])

        # legend stuff
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        jlines = pyimport("matplotlib.lines")
        man_marks = [jlines.Line2D([0], [0], color="black", linewidth=0, marker="o"),
                        jlines.Line2D([0], [0], color="black", linewidth=0, marker="D")]
        man_lbls = ["Thermal escape", "Thermal + non-thermal escape"]
        append!(man_marks, handles)
        append!(man_lbls, labels)
        ax.legend(man_marks, man_lbls, fontsize=20) 
    end

    # if moredh
    #     loss_1p5 = [h2olost(current_inventory, 1.5, 1.275, f) for f in farray]
    #     loss_10 = [h2olost(current_inventory, 10, 1.275, f) for f in farray]
    #     ax.semilogx(farray, loss_1p5, zorder=2, color=deemph)
    #     ax.semilogx(farray, loss_10, zorder=2, color=deemph)
    #     text(0.0000135, loss_1p5[2], L"(D/H)$_{M}$ = $1.5$(D/H)$_{E}$", 
    #          color=deemph, fontsize=22)
    #     text(0.0000135, loss_10[2], L"(D/H)$_{M}$ = $10$(D/H)$_{E}$", 
    #          color=deemph, fontsize=22)
    #     text(0.0000135, 90, L"(D/H)$_{M}$ = $5.5$(D/H)$_{E}$", color=emph, fontsize=22)
    # end
    
    ylim(0, 240)
    savefig("../Results/ALL STUDY PLOTS/h2oloss_$(filename_tag)m_current.png", 
            bbox_inches="tight")
    println("Water loss values for $(current_inventory):")
    println("Thermal+Nonthermal: ", loss_both)
    println("Thermal only: ", loss_thermal)
    println()
end

function plot_simple_water_loss(current_inventory=25)
    #=
    Plot just the Rayleigh distillation equation, no fractionation factor results.
    For demonstration purposes.

    current_inventory: assumed water inventory in m GEL.
    =#
    farray = range(0.000005, stop=1, length=100)

    l1 = [h2olost(current_inventory, 5.5, 1.275, f) for f in farray]

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

println("ALERT: Make sure the f values are correct. They must be filled in manually")
f_thermal = [0.0012, 0.0012, 0.016, 0.000043, 0.0022, 0.00055, 0.0013,
                 0.0011, 0.0011]
f_nonthermal = [0.038, 0.032, 0.065, 0.038, 0.068, 0.023, 0.037,
                0.037, 0.039]

println("Thermal f values: ", f_thermal)
println("Thermal+nonthermal f values: ", f_both)

plot_water_loss(f_thermal, f_nonthermal, current_inventory=20)
plot_water_loss(f_thermal, f_nonthermal, current_inventory=25)
plot_water_loss(f_thermal, f_nonthermal, current_inventory=30)

plot_water_loss(f_thermal, f_nonthermal, current_inventory=[20,25,30])