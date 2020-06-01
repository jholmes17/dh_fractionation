################################################################################
# plot_rxn_rates.jl
# TYPE: SUPPORTING (plot maker)
# WHICH: Equilibrium experiments
# DESCRIPTION: Plots chemical reaction rates by altitude and reaction.

# Eryn Cangi
# 5 April 2019
# Last edited: 3 January 2020 (NOT DONE)
# Currently tested for Julia: 0.7
################################################################################
using PyPlot
using HDF5
using LaTeXStrings
using PyCall
using PlotUtils
using JLD
using Analysis

include("PARAMETERS.jl")

threebody(k0, kinf) = :($k0*M/(1+$k0*M/$kinf)*0.6^((1+(log10($k0*M/$kinf))^2)^-1.0))
threebodyca(k0, kinf) = :($k0/(1+$k0/($kinf/M))*0.6^((1+(log10($k0/($kinf*M)))^2)^-1.0))

# NOTE: You cannot use reactionnet as defined in PARAMETERS.jl. It has to be
# this one here because we must change all the operators to be array operators.
# Net last copied: 5 May 2020
reactionnet = [   #Photodissociation
             [[:CO2], [:CO, :O], :JCO2toCOpO],
             [[:CO2], [:CO, :O1D], :JCO2toCOpO1D],
             [[:O2], [:O, :O], :JO2toOpO],
             [[:O2], [:O, :O1D], :JO2toOpO1D],
             [[:O3], [:O2, :O], :JO3toO2pO],
             [[:O3], [:O2, :O1D], :JO3toO2pO1D],
             [[:O3], [:O, :O, :O], :JO3toOpOpO],
             [[:H2], [:H, :H], :JH2toHpH],
             [[:HD], [:H, :D], :JHDtoHpD],
             [[:OH], [:O, :H], :JOHtoOpH],
             [[:OH], [:O1D, :H], :JOHtoO1DpH],
             [[:OD], [:O, :D], :JODtoOpD],
             [[:OD], [:O1D, :D], :JODtoO1DpD],
             [[:HO2], [:OH, :O], :JHO2toOHpO], # other branches should be here, but
                                               # have not been measured
             [[:DO2], [:OD, :O], :JDO2toODpO],
             [[:H2O], [:H, :OH], :JH2OtoHpOH],
             [[:H2O], [:H2, :O1D], :JH2OtoH2pO1D],
             [[:H2O], [:H, :H, :O], :JH2OtoHpHpO],
             [[:HDO], [:H, :OD], :JHDOtoHpOD], 
             [[:HDO], [:D, :OH], :JHDOtoDpOH], 
             [[:HDO], [:HD, :O1D], :JHDOtoHDpO1D], # inspiration from Yung89
             [[:HDO], [:H, :D, :O], :JHDOtoHpDpO], # inspiration from Yung89
             [[:H2O2], [:OH, :OH], :JH2O2to2OH],
             [[:H2O2], [:HO2, :H], :JH2O2toHO2pH],
             [[:H2O2], [:H2O, :O1D], :JH2O2toH2OpO1D],
             [[:HDO2], [:OH, :OD], :JHDO2toOHpOD], # Yung89
             [[:HDO2], [:DO2, :H], :JHDO2toDO2pH],
             [[:HDO2], [:HO2, :D], :JHDO2toHO2pD],
             [[:HDO2], [:HDO, :O1D], :JHDO2toHDOpO1D],

             # recombination of O
             [[:O, :O, :M], [:O2, :M], :(1.8*3.0e-33*(300 ./ T).^3.25)],
             [[:O, :O2, :N2], [:O3, :N2], :(5e-35*exp.(724 ./ T))],
             [[:O, :O2, :CO2], [:O3, :CO2], :(2.5*6.0e-34*(300 ./ T).^2.4)],
             [[:O, :O3], [:O2, :O2], :(8.0e-12*exp.(-2060 ./ T))],  # Sander 2011
             [[:O, :CO, :M], [:CO2, :M], :(2.2e-33*exp.(-1780 ./ T))],

             # O1D attack
             [[:O1D, :O2], [:O, :O2], :(3.2e-11*exp.(70 ./ T))], # verified NIST 4/3/18
             [[:O1D, :O3], [:O2, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :O3], [:O, :O, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :CO2], [:O, :CO2], :(7.5e-11*exp.(115 ./ T))], # Sander2011. NIST: 7.41e-11*exp.(120/T)
             ## O1D + H2
             [[:O1D, :H2], [:H, :OH], :(1.2e-10)],  # Sander2011. Yung89: 1e-10; NIST 1.1e-10
             [[:O1D, :HD], [:H, :OD], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1.3e-10 @298K
             [[:O1D, :HD], [:D, :OH], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1e-10 @298K
             ## O1D + H2O
             [[:O1D, :H2O], [:OH, :OH], :(1.63e-10*exp.(60 ./ T))], # Sander2011. Yung89: 2.2e-10; NIST: 1.62e-10*exp.(65/T)
             [[:O1D, :HDO], [:OD, :OH], :(1.63e-10*exp.(60 ./ T))], # Yung88: rate same as H-ana.

             # loss of H2
             [[:H2, :O], [:OH, :H], :(6.34e-12*exp.(-4000 ./ T))], # KIDA <-- Baulch, D. L. 2005
             [[:HD, :O], [:OH, :D], :(4.40e-12*exp.(-4390 ./ T))], # NIST
             [[:HD, :O], [:OD, :H], :(1.68e-12*exp.(-4400 ./ T))], # NIST
             # HD and H2 exchange
             [[:H, :HD], [:H2, :D], :(6.31e-11*exp.(-4038 ./ T))], # rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11*exp.(-3821 ./ T))], # NIST (1986, 200-300K): 8.19e-13*exp.(-2700/T)

             ## OH + H2
             [[:OH, :H2], [:H2O, :H], :(2.8e-12*exp.(-1800 ./ T))], # Sander2011. Yung89: 5.5e-12*exp.(-2000/T). KIDA: 7.7E-12*exp.(-2100/T). old rate from Mike: 9.01e-13*exp.(-1526/T)
             [[:OH, :HD], [:HDO, :H], :((3 ./ 20.)*2.8e-12*exp.(-1800 ./ T))], # Yung88: rate (3/20)*H-ana. Sander2011: 5e-12*exp.(-2130 ./ T)
             [[:OH, :HD], [:H2O, :D], :((3 ./ 20.)*2.8e-12*exp.(-1800 ./ T))], # see prev line
             [[:OD, :H2], [:HDO, :H], :(2.8e-12*exp.(-1800 ./ T))], # Yung88: rate same as H-ana (assumed)
             [[:OD, :H2], [:H2O, :D], :(0)], # Yung88 (assumed)
             ### [[:OD, :HD], [:HDO, :D], :(???)],  # possibilities for which I 
             ### [[:OD, :HD], [:D2O, :H], :(???)],  # can't find a rate...?

             # recombination of H. Use EITHER the first line OR the 2nd and 3rd.
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32*(298 ./ T).^2.27)],
             [[:H, :H, :M], [:H2, :M], :(1.6e-32*(298 ./ T).^2.27)], # general version of H+H+CO2, rate: Justin Deighan.
             [[:H, :D, :M], [:HD, :M], :(1.6e-32*(298 ./ T).^2.27)], # Yung88: rate same as H-ana.

             [[:H, :OH, :CO2], [:H2O, :CO2], :(1.9*6.8e-31*(300 ./ T).^2)], # Can't find in databases. Mike's rate.
             [[:H, :OD, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ T).^2)], # not in Yung88. assumed rate
             [[:D, :OH, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ T).^2)], # not in Yung88. assumed rate

             ## H + HO2
             [[:H, :HO2], [:OH, :OH], :(7.2e-11)], # Sander2011. Indep of T for 245<T<300
             [[:H, :HO2], [:H2, :O2], :(0.5*6.9e-12)], # 0.5 is from Krasnopolsky suggestion to Mike
             [[:H, :HO2], [:H2O, :O1D], :(1.6e-12)], # O1D is theoretically mandated
             [[:H, :DO2], [:OH, :OD], :(7.2e-11)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HD, :O2], :(0.5*6.9e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HDO, :O1D], :(1.6e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18. Yung88 has this as yielding HDO and O, not HDO and O1D
             [[:D, :HO2], [:OH, :OD], :(0.71*7.2e-11)], # Yung88: rate 0.71*H-ana (assumed). verified Yung89 3/28/18 (base: 7.05, minor disagreement)
             [[:D, :HO2], [:HD, :O2], :(0.71*0.5*6.9e-12)], # Yung88: rate 0.71*H-ana (assumed). verified Yung89 3/28/18 (base 7.29, minor disagreement)
             [[:D, :HO2], [:HDO, :O1D], :(0.71*1.6e-12)], # Yung88: rate 0.71*H-ana (assumed). Changed to O1D to match what Mike put in 3rd line from top of this section.
             [[:H, :DO2], [:HO2, :D], :(1e-10/(0.54*exp.(890 ./ T)))], # Yung88 (assumed) - turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 - turn off for Case 2

             ## H + H2O2. deuterated analogues added 3/29
             [[:H, :H2O2], [:HO2, :H2],:(2.81e-12*exp.(-1890 ./ T))], # verified NIST 4/3/18. Only valid for T>300K. No exp.eriment for lower.
             # [[:H, :HDO2], [:DO2, :H2], :(0)], # Cazaux2010: branching ratio = 0
             # [[:H, :HDO2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:DO2, :H2], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             [[:H, :H2O2], [:H2O, :OH],:(1.7e-11*exp.(-1800 ./ T))], # verified NIST 4/3/18
             [[:H, :HDO2], [:HDO, :OH], :(0.5*1.16e-11*exp.(-2110 ./ T))], # Cazaux2010: BR = 0.5. Rate for D + H2O2, valid 294<T<464K, NIST, 4/3/18
             [[:H, :HDO2], [:H2O, :OD], :(0.5*1.16e-11*exp.(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:HDO, :OH], :(0.5*1.16e-11*exp.(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:H2O, :OD], :(0.5*1.16e-11*exp.(-2110 ./ T))], # see previous line
             [[:D, :HDO2], [:OD, :HDO], :(0.5*1.16e-11*exp.(-2110 ./ T))], # added 4/3 with assumed rate from other rxns
             [[:D, :HDO2], [:OH, :D2O], :(0.5*1.16e-11*exp.(-2110/T))], # sourced from Cazaux et al

             # Interconversion of odd H
             ## H + O2
             [[:H, :O2], [:HO2], threebody(:(2.0*4.4e-32*(T/300.).^-1.3), # Sander2011, 300K+. Yung89: 5.5e-32(T/300).^-1.6, 7.5e-11 valid 200-300K.
                                           :(7.5e-11*(T/300.).^0.2))],  # NIST has the temp info.
             [[:D, :O2], [:DO2], threebody(:(2.0*4.4e-32*(T/300.).^-1.3), # Yung88: rate same as H-ana.
                                           :(7.5e-11*(T/300.).^0.2))],

             ## H + O3
             [[:H, :O3], [:OH, :O2], :(1.4e-10*exp.(-470 ./ T))], # verified Yung89, NIST 4/3/18
             [[:D, :O3], [:OD, :O2], :(0.71*1.4e-10*exp.(-470 ./ T))], # Yung88: rate 0.71*H-ana (assumed). verified Yung89, NIST 4/3/18.
             ## O + OH
             [[:O, :OH], [:O2, :H], :(1.8e-11*exp.(180 ./ T))], # Sander2011. KIDA+NIST 4/3/18 150-500K: 2.4e-11*exp.(110 ./ T). Yung89: 2.2e-11*exp.(120/T) for both this and D analogue.
             [[:O, :OD], [:O2, :D], :(1.8e-11*exp.(180 ./ T))], # Yung88: rate same as H-ana.
             ## O + HO2
             [[:O, :HO2], [:OH, :O2], :(3.0e-11*exp.(200 ./ T))], # Sander2011. KIDA (220-400K): 2.7e-11*exp.(224/T)
             [[:O, :DO2], [:OD, :O2], :(3.0e-11*exp.(200 ./ T))], # Yung88: rate same as H-ana. verified Yung89 4/3/18
             ## O + H2O2
             [[:O, :H2O2], [:OH, :HO2], :(1.4e-12*exp.(-2000 ./ T))], # Sander2011. verified NIST 4/3/18.
             [[:O, :HDO2], [:OD, :HO2], :(0.5*1.4e-12*exp.(-2000 ./ T))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             [[:O, :HDO2], [:OH, :DO2], :(0.5*1.4e-12*exp.(-2000 ./ T))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             ## OH + OH
             [[:OH, :OH], [:H2O, :O], :(4.2e-12*exp.(-240 ./ T))], # NIST+KIDA, 200-350K: 6.2e-14*(T/300).^2.62*exp.(945 ./ T) changed 4/3/18. Yung89: 4.2e-12*exp.(-240/T). old rate w/mystery origin: 1.8e-12.
             [[:OD, :OH], [:HDO, :O], :(4.2e-12*exp.(-240 ./ T))], # Yung88: rate same as H-ana
             [[:OH, :OH], [:H2O2], threebody(:(1.3*6.9e-31*(T/300.).^-1.0),:(2.6e-11))], # Sander2011. Why 1.3?
             [[:OD, :OH], [:HDO2], threebody(:(1.3*6.9e-31*(T/300.).^-1.0),:(2.6e-11))], # Yung88: rate same as H-ana
             ## OH + O3
             [[:OH, :O3], [:HO2, :O2], :(1.7e-12*exp.(-940 ./ T))], # Sander2011, temp by NIST 220-450K. Yung89: 1.6 not 1.7 -> temp 200-300K by NIST (older info)
             [[:OD, :O3], [:DO2, :O2], :(1.7e-12*exp.(-940 ./ T))], # Yung88: rate same as H-ana
             ## OH + HO2
             [[:OH, :HO2], [:H2O, :O2], :(4.8e-11*exp.(250 ./ T))], # verified NIST 4/3/18. Yung89: 4.6e-11*exp.(230/T) for this and next 2.
             [[:OH, :DO2], [:HDO, :O2], :(4.8e-11*exp.(250 ./ T))], # Yung88: same as H-ana.
             [[:OD, :HO2], [:HDO, :O2], :(4.8e-11*exp.(250 ./ T))], # Yung88: same as H-ana.
             ## OH + H2O2
             [[:OH, :H2O2], [:H2O, :HO2], :(2.9e-12*exp.(-160 ./ T))], # NIST+KIDA 4/3/18, valid 240-460K. Yung89: 3.3e-12*exp.(-200/T). Sander2011 recommends an average value of 1.8e-12, but this seems too high for martian temps
             [[:OD, :H2O2], [:HDO, :HO2], :(2.9e-12*exp.(-160 ./ T))], # Yung88: same as H-ana (assumed)
             [[:OD, :H2O2], [:H2O, :DO2], :(0)],  # Yung88 (assumed)
             [[:OH, :HDO2], [:HDO, :HO2], :(0.5*2.9e-12*exp.(-160 ./ T))], # Yung88: rate 0.5*H-ana.
             [[:OH, :HDO2], [:H2O, :DO2], :(0.5*2.9e-12*exp.(-160 ./ T))], # Yung88: rate 0.5*H-ana.
             ## HO2 + O3
             [[:HO2, :O3], [:OH, :O2, :O2], :(1.0e-14*exp.(-490 ./ T))], # Sander2011. Yung89: 1.1e-14*exp.(-500/T). KIDA 250-340K: 2.03e-16*(T/300).^4.57*exp.(693/T). All give comparable rate values (8.6e-16 to 1e-15 at 200K)
             [[:DO2, :O3], [:OD, :O2, :O2], :(1.0e-14*exp.(-490 ./ T))], # Yung88: same as H-ana (assumed)
             ## HO2 + HO2
             [[:HO2, :HO2], [:H2O2, :O2], :(3.0e-13*exp.(460 ./ T))], # Sander2011. Yung89: 2.3e-13*exp.(600/T). KIDA 230-420K: 2.2e-13*exp.(600/T)
             [[:DO2, :HO2], [:HDO2, :O2], :(3.0e-13*exp.(460 ./ T))], # Yung88: same as H-ana (assumed)
             # *** why do we have he next two reactions? I forgot...
             [[:HO2, :HO2, :M], [:H2O2, :O2, :M], :(2*2.1e-33*exp.(920 ./ T))], # Sander2011.
             [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :(2*2.1e-33*exp.(920 ./ T))], # added 3/13 with assumed same rate as H analogue

             ## OH + D or OD + H (no non-deuterated analogues)
             [[:OD, :H], [:OH, :D], :(3.3e-9*(T.^-0.63)/(0.72*exp.(717 ./ T)))], # rate: Yung88. NIST (Howard82): 5.25E-11*(T/298).^-0.63  - turn off for Case 2
             [[:OH, :D], [:OD, :H], :(3.3e-9*T.^-0.63)], # Yung88  - turn off for Case 2

             # CO2 recombination due to odd H (with HOCO intermediate)
             ## straight to CO2
             [[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13*(T/300.).^0.6),:(2.1e9*(T/300.).^6.1))], # Sander2011
             [[:CO, :OD], [:CO2, :D], threebodyca(:(1.5e-13*(T/300.).^0.6),:(2.1e9*(T/300.).^6.1))], # Yung88: same as H-ana.
             ### possible deuterated analogues below
             [[:OH, :CO], [:HOCO], threebody(:(5.9e-33*(T/300.).^-1.4),:(1.1e-12*(T/300.).^1.3))], # Sander2011
             [[:OD, :CO], [:DOCO], threebody(:(5.9e-33*(T/300.).^-1.4),:(1.1e-12*(T/300.).^1.3))],

             [[:HOCO, :O2], [:HO2, :CO2], :(2.09e-12)], # verified NIST 4/3/18
             [[:DOCO, :O2], [:DO2,:CO2], :(2.09e-12)],  # assumed?

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 / Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5)*8.7e-10)]
             ]

function make_ratexdensity(n_current, temparr, species, species_role)
    #=
    n_current: a given result file for a converged atmosphere
    temparr: Array of temperatures by altitude for the atmosphere
    species: only reactions including this species will be plotted
    species_role: whether to lo look for the species as a reactant, product, or both
    =#

    rxn_dat =  Dict{String,Array{Float64, 1}}()

    for rxn in reactionnet
        reactants = rxn[1]
        products = rxn[2]

        role = Dict("reactant"=>reactants, "product"=>products, 
                    "both"=>[reactants, products])

        if ~in(species, role[species_role])
            continue
        end

        # get the reactants and products in string form for use in plot labels
        reacts = join(rxn[1], " + ")
        prods = join(rxn[2], " + ")
        rxn_str = string(reacts) * " --> " * string(prods)

        # calculate the reaction strength, which is the rate x density of least 
        # dense species at given altitude.
        if typeof(rxn[3]) == Symbol # for photodissociation
            alt_arr = n_current[rxn[1][1]]
            rate_arr = n_current[rxn[3]]
            ratexdens = rate_arr .* alt_arr
            # put ratexdens in dictionary - photodissociation has only 1 reactant
            rxn_dat[rxn_str] = ratexdens
        else  
            # for reactions with more than one reactant. the procedure is to 
            # multiply
            
            lowest_ratexdens = fill(9e50, (length(temparr)))
            for r in rxn[1]
                if r != :M
                    # species densities by altitude 
                    alt_arr = n_current[r]

                    # reaction rate: this is witchcraft and I forgot how I figured it out.
                    rate = rxn[3]           # for regular rxns
                    @eval ratefunc(T) = $rate
                    rate_arr = Base.invokelatest(ratefunc, temparr)
                    ratexdens = rate_arr .* alt_arr
                    
                    # reactions are limited by the lowest-density value, so must look for lowest value of ratexdens across reactants.
                    # i verified that this works by looking at the ratexdens for O and O2 at element 39 in rxn O + O2 + N2.
                    for i in range(1, length=length(ratexdens))
                        if ratexdens[i] < lowest_ratexdens[i]
                            lowest_ratexdens[i] = ratexdens[i]
                        end
                    end

                end
            end

            # set all values of 9e50 to 0 since there was no lowest
            for iter in eachindex(lowest_ratexdens)
                if lowest_ratexdens[iter] == 9e50
                    lowest_ratexdens[iter] = 0
                end
            end

            rxn_dat[rxn_str] = lowest_ratexdens
        end
    end
    return rxn_dat
end

# ==============================================================================

basepath = results_dir*"TradeoffPlots/Main Results/"

# Do stuff for surface temperature 
for Ts in collect(150:10:270)

    # set up readfile and temperature array
    readfile = basepath*"temp_$(Ts)_$(meanTtint)_$(meanTeint)/converged_temp_$(Ts)_$(meanTtint)_$(meanTeint).h5"
    global alt=h5read(readfile,"n_current/alt")

    # make the temperature array
    tempsbyalt = Array{Float64}(undef, length(alt)-2)

    Temp(z::Float64) = Tpiecewise(z, Ts, meanTt, meanTe, "surf")
    i = 1
    for i in range(1, length=length(alt)-2)
        tempsbyalt[i] = Temp(alt[i])
        i += 1
    end

    # retrieve the matrix and make the rate x density plot
    ncur = get_ncurrent(readfile)
    species_of_interest = :O3
    role = "product"

    rxd = make_ratexdensity(ncur, tempsbyalt, species_of_interest, role)

    THRESHHOLD = 1e-50  # A given reaction must have a rate above this value at 
                       # some point in the atmosphere to be plotted
    fig, ax = subplots(figsize=(9,7))
    subplots_adjust(wspace=0, bottom=0.25)
    plot_bg(ax)
    rxnlist = []
    surfvallist = []
    for kv in rxd
        lbl = "$(kv[1])"
        if any(x->x>THRESHHOLD, kv[2][1]) # this only counts reactions with a rate over 1.
            push!(rxnlist, kv[1])
            append!(surfvallist, kv[2][1])
        end
    end
    bar(rxnlist, surfvallist, zorder=10)
    title("Reaction rates at the surface")
    xlabel("Reaction")
    xticks(rotation=45)
    yscale("log")
    ylim(1, 1e8)
    ylabel("Reaction rate (#/cm^3/s) (check this)")
    savefig(basepath*"reaction plots"*"/rxns_surface_$(Ts).png", bbox_inches="tight")
    close(fig)


    fig2, ax2 = subplots(figsize=(9,7))
    subplots_adjust(wspace=0, bottom=0.25)
    plot_bg(ax2)
    for kv in rxd
        lw = 1
        ls = "-"
        lbl = "$(kv[1])"
        if any(x->x>THRESHHOLD, kv[2])
            semilogx(kv[2], alt[1:length(alt)-2]./1e5, linestyle=ls, linewidth=lw, label=lbl)
        end
    end
    title("Reaction rates vs. altitude")
    # xlim(1e-12, 1e8)
    ylabel("Altitude (km)")
    xlabel("Reaction rate (#/cm^3/s) (check this)")
    legend()
    savefig(basepath*"reaction plots/"*"$(species_of_interest)_$(role)_reactions_Ts$(Ts).png", bbox_inches="tight")
    close(fig)
end

