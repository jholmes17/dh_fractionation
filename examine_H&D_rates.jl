################################################################################
# examine_H&D_rates.jl
# TYPE: Supporting (Troubleshooting)
# WHICH: Equilibrium and perturbation experiments
# DESCRIPTION: find rate profiles by altitude for H and D production.
#
# Eryn Cangi
# 3 June 2019
# Currently tested for Julia: 0.7
################################################################################

using PyPlot
using HDF5
using LaTeXStrings
using PyCall
using PlotUtils
using JLD

# Constants and reaction net ===================================================
alt = 0:2e5:200e5
const intaltgrid = round.(Int64, alt/1e5)[2:end-1];
threebody(k0, kinf) = :($k0*M/(1+$k0*M/$kinf)*0.6^((1+(log10($k0*M/$kinf))^2)^-1))
threebodyca(k0, kinf) = :($k0/(1+$k0/($kinf/M))*0.6^((1+(log10($k0/($kinf*M)))^2)^-1))

reactionnet = [   #Photodissociation
             [[:CO2], [:CO, :O], :JCO2toCOpO],
             [[:CO2], [:CO, :O1D], :JCO2toCOpO1D],
             [[:O2], [:O, :O], :JO2toOpO],
             [[:O2], [:O, :O1D], :JO2toOpO1D],
             [[:O3], [:O2, :O], :JO3toO2pO],
             [[:O3], [:O2, :O1D], :JO3toO2pO1D],
             [[:O3], [:O, :O, :O], :JO3toOpOpO],
             [[:H2], [:H, :H], :JH2toHpH],
             [[:HD], [:H, :D], :JHDtoHpD],  # added 3/29 using H2 xsects
             [[:OH], [:O, :H], :JOHtoOpH],
             [[:OH], [:O1D, :H], :JOHtoO1DpH],
             [[:OD], [:O, :D], :JODtoOpD], # added 3/28
             [[:OD], [:O1D, :D], :JODtoO1DpD],
             [[:HO2], [:OH, :O], :JHO2toOHpO], # other branches should be here, but
                                               # have not been measured
             [[:DO2], [:OD, :O], :JDO2toODpO],  # added 3/29 using HO2 xsects
             [[:H2O], [:H, :OH], :JH2OtoHpOH],
             [[:H2O], [:H2, :O1D], :JH2OtoH2pO1D],
             [[:H2O], [:H, :H, :O], :JH2OtoHpHpO],
             [[:HDO], [:H, :OD], :JHDOtoHpOD], 
             [[:HDO], [:D, :OH], :JHDOtoDpOH], 
             [[:HDO], [:HD, :O1D], :JHDOtoHDpO1D], # added 3/28, inspiration from Yung89
             [[:HDO], [:H, :D, :O], :JHDOtoHpDpO], # added 3/28, inspiration from Yung89
             [[:H2O2], [:OH, :OH], :JH2O2to2OH],
             [[:H2O2], [:HO2, :H], :JH2O2toHO2pH],
             [[:H2O2], [:H2O, :O1D], :JH2O2toH2OpO1D],
             [[:HDO2], [:OH, :OD], :JHDO2toOHpOD], #  Yung89; using xsect for H2O2
             [[:HDO2], [:DO2, :H], :JHDO2toDO2pH], # added 3/30 using H2O2 xsects
             [[:HDO2], [:HO2, :D], :JHDO2toHO2pD], # added 3/30 using H2O2 xsects
             [[:HDO2], [:HDO, :O1D], :JHDO2toHDOpO1D], # added 3/30 using H2O2 xsects

             # recombination of O
             [[:O, :O, :M], [:O2, :M], :(1.8*3.0e-33*(300 ./ T)^3.25)],
             [[:O, :O2, :N2], [:O3, :N2], :(5e-35*exp(724 ./ T))],
             [[:O, :O2, :CO2], [:O3, :CO2], :(2.5*6.0e-34*(300 ./ T)^2.4)],
             [[:O, :O3], [:O2, :O2], :(8.0e-12*exp(-2060 ./ T))],
             [[:O, :CO, :M], [:CO2, :M], :(2.2e-33*exp(-1780 ./ T))],

             # O1D attack
             [[:O1D, :O2], [:O, :O2], :(3.2e-11*exp(70 ./ T))], # verified NIST 4/3/18
             [[:O1D, :O3], [:O2, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :O3], [:O, :O, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :CO2], [:O, :CO2], :(7.5e-11*exp(115 ./ T))], # Sander2011. NIST: 7.41e-11*exp(120/T)
             ## O1D + H2
             [[:O1D, :H2], [:H, :OH], :(1.2e-10)],  # Sander2011. Yung89: 1e-10; NIST 1.1e-10
             [[:O1D, :HD], [:H, :OD], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1.3e-10 @298K
             [[:O1D, :HD], [:D, :OH], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1e-10 @298K
             ## O1D + H2O
             [[:O1D, :H2O], [:OH, :OH], :(1.63e-10*exp(60 ./ T))], # Sander2011. Yung89: 2.2e-10; NIST: 1.62e-10*exp(65/T)
             [[:O1D, :HDO], [:OD, :OH], :(1.63e-10*exp(60 ./ T))], # Yung88: rate  same as H-ana.

             # loss of H2
             [[:H2, :O], [:OH, :H], :(6.34e-12*exp(-4000 ./ T))], # Mike's rate? I believe this based on next 2 lines. ~10^-20 at 200K
             [[:HD, :O], [:OH, :D], :(4.40e-12*exp(-4390 ./ T))], # not in Yung. rate: NIST. ~10^-21 at 200K
             [[:HD, :O], [:OD, :H], :(1.68e-12*exp(-4400 ./ T))], # not in Yung. rate: NIST. ~10^-22 at 200K
             # HD and H2 exchange
             [[:H, :HD], [:H2, :D], :(6.31e-11*exp(-4038 ./ T))], # rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11*exp(-3821 ./ T))], # NIST (1986, 200-300K): 8.19e-13*exp(-2700/T)

             ## OH + H2
             [[:OH, :H2], [:H2O, :H], :(2.8e-12*exp(-1800 ./ T))], # Sander2011. Yung89: 5.5e-12*exp(-2000/T). KIDA: 7.7E-12*exp(-2100/T). old rate from Mike: 9.01e-13*exp(-1526/T)
             [[:OH, :HD], [:HDO, :H], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ T))], # Yung88: rate (3/20)*H-ana. Sander2011: 5e-12*exp(-2130 ./ T)
             [[:OH, :HD], [:H2O, :D], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ T))], # see prev line
             [[:OD, :H2], [:HDO, :H], :(2.8e-12*exp(-1800 ./ T))], # Yung88: rate same as H-ana (assumed)
             [[:OD, :H2], [:H2O, :D], :(0)], # Yung88 (assumed)
             ### [[:OD, :HD], [:HDO, :D], :(???)],
             ### [[:OD, :HD], [:D2O, :H], :(???)],

             # recombination of H
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32*(298 ./ T)^2.27)], # keep this off unless next 2 off
             [[:H, :H, :M], [:H2, :M], :(1.6e-32*(298 ./ T)^2.27)], # general version of H+H+CO2, rate: Justin Deighan.
             [[:H, :D, :M], [:HD, :M], :(1.6e-32*(298 ./ T)^2.27)], # Yung88: rate same as H-ana.

             [[:H, :OH, :CO2], [:H2O, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # Can't find in databases. Mike's rate.
             [[:H, :OD, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # not in Yung88. assumed rate
             [[:D, :OH, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # not in Yung88. assumed rate

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
             [[:H, :DO2], [:HO2, :D], :(1e-10/(0.54*exp(890 ./ T)))], # Yung88  (assumed) - turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 - turn off for Case 2

             ## H + H2O2. deuterated analogues added 3/29
             [[:H, :H2O2], [:HO2, :H2],:(2.81e-12*exp(-1890 ./ T))], # verified NIST 4/3/18. Only valid for T>300K. No experiment for lower.
             [[:H, :HDO2], [:DO2,:H2], :(0)], # Cazaux2010: BR = 0
             [[:H, :HDO2], [:HO2,:HD], :(0)], # Cazaux2010: BR = 0
             [[:D, :H2O2], [:DO2,:H2], :(0)], # Cazaux2010: BR = 0
             [[:D, :H2O2], [:HO2,:HD], :(0)], # Cazaux2010: BR = 0
             [[:H, :H2O2], [:H2O, :OH],:(1.7e-11*exp(-1800 ./ T))], # verified NIST 4/3/18
             [[:H, :HDO2], [:HDO,:OH], :(0.5*1.16E-11*exp(-2110 ./ T))], # Cazaux2010: BR = 0.5. Rate for D + H2O2, valid 294<T<464K, NIST, 4/3/18
             [[:H, :HDO2], [:H2O,:OD], :(0.5*1.16E-11*exp(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:HDO,:OH], :(0.5*1.16E-11*exp(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:H2O,:OD], :(0.5*1.16E-11*exp(-2110 ./ T))], # see previous line
             [[:D, :HDO2], [:OD, :HDO], :(0.5*1.16E-11*exp(-2110 ./ T))], # added 4/3 with assumed rate from other rxns
             [[:D, :HDO2], [:OH, :D2O], :(0.5*1.16E-11*exp(-2110/T))], # sourced from Cazaux et al

             # Interconversion of odd H
             ## H + O2
             [[:H, :O2], [:HO2], threebody(:(2.0*4.4e-32*(T/300.)^-1.3), # Sander2011, 300K+. Yung89: 5.5e-32(T/300)^-1.6, 7.5e-11 valid 200-300K.
                                           :(7.5e-11*(T/300.)^0.2))],  # NIST has the temp info.
             [[:D, :O2], [:DO2], threebody(:(2.0*4.4e-32*(T/300.)^-1.3), # Yung88: rate same as H-ana.
                                           :(7.5e-11*(T/300.)^0.2))],

             ## H + O3
             [[:H, :O3], [:OH, :O2], :(1.4e-10*exp(-470 ./ T))], # verified Yung89, NIST 4/3/18
             [[:D, :O3], [:OD, :O2], :(0.71*1.4e-10*exp(-470 ./ T))], # Yung88: rate 0.71*H-ana (assumed). verified Yung89, NIST 4/3/18.
             ## O + OH
             [[:O, :OH], [:O2, :H], :(1.8e-11*exp(180 ./ T))], # Sander2011. KIDA+NIST 4/3/18 150-500K: 2.4e-11*exp(110 ./ T). Yung89: 2.2e-11*exp(120/T) for both this and D analogue.
             [[:O, :OD], [:O2, :D], :(1.8e-11*exp(180 ./ T))], # Yung88: rate same as H-ana.
             ## O + HO2
             [[:O, :HO2], [:OH, :O2], :(3.0e-11*exp(200 ./ T))], # Sander2011. KIDA (220-400K): 2.7e-11*exp(224/T)
             [[:O, :DO2], [:OD, :O2], :(3.0e-11*exp(200 ./ T))], # Yung88: rate same as H-ana. verified Yung89 4/3/18
             ## O + H2O2
             [[:O, :H2O2], [:OH, :HO2], :(1.4e-12*exp(-2000 ./ T))], # Sander2011. verified NIST 4/3/18.
             [[:O, :HDO2], [:OD, :HO2], :(0.5*1.4e-12*exp(-2000 ./ T))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             [[:O, :HDO2], [:OH, :DO2], :(0.5*1.4e-12*exp(-2000 ./ T))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             ## OH + OH
             [[:OH, :OH], [:H2O, :O], :(4.2e-12*exp(-240 ./ T))], # NIST+KIDA, 200-350K: 6.2e-14*(T/300)^2.62*exp(945 ./ T) changed 4/3/18. Yung89: 4.2e-12*exp(-240/T). old rate w/mystery origin: 1.8e-12.
             [[:OD, :OH], [:HDO, :O], :(4.2e-12*exp(-240 ./ T))], # Yung88: rate same as H-ana
             [[:OH, :OH], [:H2O2], threebody(:(1.3*6.9e-31*(T/300.)^-1.0),:(2.6e-11))], # Sander2011. Why 1.3?
             [[:OD, :OH], [:HDO2], threebody(:(1.3*6.9e-31*(T/300.)^-1.0),:(2.6e-11))], # Yung88: rate same as H-ana
             ## OH + O3
             [[:OH, :O3], [:HO2, :O2], :(1.7e-12*exp(-940 ./ T))], # Sander2011, temp by NIST 220-450K. Yung89: 1.6 not 1.7 -> temp 200-300K by NIST (older info)
             [[:OD, :O3], [:DO2, :O2], :(1.7e-12*exp(-940 ./ T))], # Yung88: rate same as H-ana
             ## OH + HO2
             [[:OH, :HO2], [:H2O, :O2], :(4.8e-11*exp(250 ./ T))], # verified NIST 4/3/18. Yung89: 4.6e-11*exp(230/T) for this and next 2.
             [[:OH, :DO2], [:HDO, :O2], :(4.8e-11*exp(250 ./ T))], # Yung88: same as H-ana.
             [[:OD, :HO2], [:HDO, :O2], :(4.8e-11*exp(250 ./ T))], # Yung88: same as H-ana.
             ## OH + H2O2
             [[:OH, :H2O2], [:H2O, :HO2], :(2.9e-12*exp(-160 ./ T))], # NIST+KIDA 4/3/18, valid 240-460K. Yung89: 3.3e-12*exp(-200/T). Sander2011 recommends an average value of 1.8e-12, but this seems too high for martian temps
             [[:OD, :H2O2], [:HDO, :HO2], :(2.9e-12*exp(-160 ./ T))], # Yung88: same as H-ana (assumed)
             [[:OD, :H2O2], [:H2O, :DO2], :(0)],  # Yung88 (assumed)
             [[:OH, :HDO2], [:HDO, :HO2], :(0.5*2.9e-12*exp(-160 ./ T))], # Yung88: rate 0.5*H-ana.
             [[:OH, :HDO2], [:H2O, :DO2], :(0.5*2.9e-12*exp(-160 ./ T))], # Yung88: rate 0.5*H-ana.
             ## HO2 + O3
             [[:HO2, :O3], [:OH, :O2, :O2], :(1.0e-14*exp(-490 ./ T))], # Sander2011. Yung89: 1.1e-14*exp(-500/T). KIDA 250-340K: 2.03e-16*(T/300)^4.57*exp(693/T). All give comparable rate values (8.6e-16 to 1e-15 at 200K)
             [[:DO2, :O3], [:OD, :O2, :O2], :(1.0e-14*exp(-490 ./ T))], # Yung88: same as H-ana (assumed)
             ## HO2 + HO2
             [[:HO2, :HO2], [:H2O2, :O2], :(3.0e-13*exp(460 ./ T))], # Sander2011. Yung89: 2.3e-13*exp(600/T). KIDA 230-420K: 2.2e-13*exp(600/T)
             [[:DO2, :HO2], [:HDO2, :O2], :(3.0e-13*exp(460 ./ T))], # Yung88: same as H-ana (assumed)
             # *** why do we have he next two reactions? I forgot...
             [[:HO2, :HO2, :M], [:H2O2, :O2, :M], :(2*2.1e-33*exp(920 ./ T))], # Sander2011.
             [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :(2*2.1e-33*exp(920 ./ T))], # added 3/13 with assumed same rate as H analogue

             ## OH + D or OD + H (no non-deuterated analogues)
             [[:OD, :H], [:OH, :D], :(3.3e-9*(T^-0.63)/(0.72*exp(717 ./ T)))], # rate: Yung88. NIST (Howard82): 5.25E-11*(T/298)^-0.63  - turn off for Case 2
             [[:OH, :D], [:OD, :H], :(3.3e-9*T^-0.63)], # Yung88  - turn off for Case 2

             # CO2 recombination due to odd H (with HOCO intermediate)
             ## straight to CO2
             [[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13*(T/300.)^0.6),:(2.1e9*(T/300.)^6.1))], # Sander2011
             [[:CO, :OD], [:CO2, :D], threebodyca(:(1.5e-13*(T/300.)^0.6),:(2.1e9*(T/300.)^6.1))], # Yung88: same as H-ana.
             ### possible deuterated analogues below
             [[:OH, :CO], [:HOCO], threebody(:(5.9e-33*(T/300.)^-1.4),:(1.1e-12*(T/300.)^1.3))], # Sander2011
             [[:OD, :CO], [:DOCO], threebody(:(5.9e-33*(T/300.)^-1.4),:(1.1e-12*(T/300.)^1.3))],

             [[:HOCO, :O2], [:HO2, :CO2], :(2.09e-12)], # verified NIST 4/3/18
             [[:DOCO, :O2], [:DO2,:CO2], :(2.09e-12)],  # assumed?

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 / Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5)*8.7e-10)],
             ];

const Jratelist=[:JCO2ion,:JCO2toCOpO,:JCO2toCOpO1D,:JO2toOpO,:JO2toOpO1D,
                 :JO3toO2pO,:JO3toO2pO1D,:JO3toOpOpO,:JH2toHpH,:JOHtoOpH,
                 :JOHtoO1DpH,:JHO2toOHpO,:JH2OtoHpOH,:JH2OtoH2pO1D,:JH2OtoHpHpO,
                 :JH2O2to2OH,:JH2O2toHO2pH,:JH2O2toH2OpO1D,
                 # deuterated species J rates:
                 :JHDOtoHpOD, :JHDOtoDpOH, :JHDO2toOHpOD,
                 # new March 2018
                 :JHDOtoHDpO1D, :JHDOtoHpDpO, :JODtoOpD, :JHDtoHpD, :JDO2toODpO,
                 :JHDO2toDO2pH, :JHDO2toHO2pD, :JHDO2toHDOpO1D, :JODtoO1DpD
                 ];

specieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
                         :O1D, :H, :N2, :Ar, :CO2pl, :HOCO,
                         # species for deuterium chemistry:
                         :HDO, :OD, :HDO2, :D, :DO2, :HD, :DOCO];

# Functions ====================================================================
@eval begin
    function reactionrates_local($(specieslist...), $(Jratelist...), T, M)
        #= a function to return chemical reaction rates for specified species
           concentrations =#
        $(Expr(:vcat, map(x->Expr(:call,:*,x[1]..., x[3]), reactionnet)...))
    end
end

function reactionrates(n_current)
    #=
    Creates an array of size length(intaltgrid) x (number of reactions).
    Populated with chemical reaction rates for each reaction based on species
    populations.
    =#
    theserates = fill(convert(Float64, NaN),(length(intaltgrid),length(reactionnet)))
    for ialt in 1:length(intaltgrid)
        theserates[ialt,:] = reactionrates_local([[n_current[sp][ialt] for sp in specieslist];
                                                [n_current[J][ialt] for J in Jratelist];
                                                Temp(alt[ialt+1]);
                                                n_tot(n_current, alt[ialt+1])]...)
    end
    rxnlist = ["$(join(x[1], " + ")) -> $(join(x[2], " + "))" for x in reactionnet]
    return theserates, rxnlist
end

function get_colors(L)
    #=
    Generates some colors based on a color map for use in plotting a bunch of
    lines all at once.
    L: number of colors to generate.
    =#

    colors_dumb = [cgrad(:viridis)[x] for x in range(0, stop=1, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = red(colors_dumb[i])
        c[i, 2] = green(colors_dumb[i])
        c[i, 3] = blue(colors_dumb[i])
    end
    return c
end

function get_ncurrent(readfile)
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

function Tpiecewise(z::Float64, Tsurf, Ttropo, Tinf, lapserate=-1.4e-5, ztropo=90e5)
    #=
    DO NOT MODIFY! If you want to change the temperature, define a
    new function or select different arguments and pass to Temp(z)

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 temperatures for altitudes
    >htropo=90km, fixed at Ttropo=125K between htropo and
    htropo-htropowidth=60km, and rising at a constant lapse rate
    (1.4K/km) below.

    z: an altitude in cm.
    returns: temperature in K at the requested altitude.
    =#
    # allow varying troposphere width
    ztropowidth = ztropo - (1/lapserate)*(Ttropo-Tsurf)
    if ztropowidth < 0   # in case a profile is nonphysical, let lapse rate vary
        ztropowidth = 30e5
        m = (ztropo/1e5 - ztropowidth/1e5) / (Ttropo - Tsurf)
        lapserate = (1/m)*1e-5
        print(lapserate)
    end

    if z >= ztropo
        return Tinf - (Tinf - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Tinf))
    end
    if ztropo > z >= ztropo - ztropowidth
        return Ttropo
    end
    if ztropo-ztropowidth > z
        return Ttropo-lapserate*(ztropo-ztropowidth-z)
    end
end
Temp(z::Float64) = Tpiecewise(z, 192.0, 110.0, 199.0)

function n_tot(n_current, z)
    # get the total number density at a given altitude
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])


# Main routine =================================================================
base = "/home/emc/GDrive-CU/Research/Results/TradeoffPlots/Oflux_1.2e8/"
ncur_dh_bulge = get_ncurrent(base*"converged_standardwater_D_Oflux_1.2e8.h5")
ratetbl, rxnlist = reactionrates(ncur_dh_bulge)
println(rxnlist)
# Find the reactions associated with H or D loss and production.
Hloss = []
Hprod = []
Dloss = []
Dprod = []

for R in rxnlist
    if occursin(r"\bH\b.+(?=->)", R)
        push!(Hloss, R)
    end
    if occursin(r"(?<=->).+\bH\b", R)
        push!(Hprod, R)
    end
    if occursin(r"\bD\b.+(?=->)", R)
        push!(Dloss, R)
    end
    if occursin(r"(?<=->).+\bD\b", R)
        push!(Dprod, R)
    end
end

# Find the indices of the loss and production reactions within rxnlist. These indices correspond to the column number in ratetbl.
Hloss_i = []
Hprod_i = []
Dloss_i = []
Dprod_i = []
for rxn in Hloss
    append!(Hloss_i, findfirst(isequal(rxn), rxnlist))
end
for rxn in Hprod
    append!(Hprod_i, findfirst(isequal(rxn), rxnlist))
end
for rxn in Dloss
    append!(Dloss_i, findfirst(isequal(rxn), rxnlist))
end
for rxn in Dprod
    append!(Dprod_i, findfirst(isequal(rxn), rxnlist))
end

Hloss_profile = Array{Float64}(undef, 99)
Hprod_profile = Array{Float64}(undef, 99)
Dloss_profile = Array{Float64}(undef, 99)
Dprod_profile = Array{Float64}(undef, 99)

for i in Hloss_i
    Hloss_profile .+= ratetbl[:, i]
end

for i in Hprod_i
    Hprod_profile .+= ratetbl[:, i]
end

for i in Dloss_i
    Dloss_profile .+= ratetbl[:, i]
end

for i in Dprod_i
    Dprod_profile .+= ratetbl[:, i]
end

Hloss_profile .*= 2e5
Hprod_profile .*= 2e5
Dloss_profile .*= 2e5
Dprod_profile .*= 2e5

# Plotting =====================================================================

marks = [".","o","v","^","<",">","8","s","*","x","X","D", ".","o","v","^",
		 "<",">","8","s","*","x","X","D"]

 # Set plot fonts to be good for posters and presentations
rcParams = PyCall.PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 22
rcParams["axes.labelsize"]= 24
rcParams["xtick.labelsize"] = 22
rcParams["ytick.labelsize"] = 22

# Total H and D loss and production rates --------------------------------------
fig, ax = subplots(figsize=(10,5))
plot(Hprod_profile, intaltgrid, color="red", label="H prod")
plot(Hloss_profile, intaltgrid, color="navy", label="H loss")
plot(Dprod_profile, intaltgrid, color="red", linestyle="--", label="D prod")
plot(Dloss_profile, intaltgrid, color="navy", linestyle="--", label="D loss")
xscale("log")
ylabel("Altitude (km)")
xlabel(L"Reaction rate (cm$^{-2}$)")
legend()
savefig(base*"total_rates_H_and_D.png", bbox_inches="tight")

# Now see which reaction seems to be responsible. H loss and production:
fig = figure(figsize=(8,6))
cullers = get_colors(length(Hloss_i))
c = 1
for i in Hloss_i
    plot(ratetbl[:, i].*2e5, intaltgrid, label=rxnlist[i], color=cullers[c, :], marker=marks[c], ms=7)
    c+=1
end
xscale("log")
legend(bbox_to_anchor=(1.1, 1))
title("Chemical reactions: H loss")
xlabel(L"Reaction rate (cm$^{-2}$)")
ylabel("Altitude (km)")
ylim(25, 75)
savefig(base*"Hloss.png", bbox_inches="tight")

fig = figure(figsize=(8,6))
cullers = get_colors(length(Hprod_i))
c = 1
for i in Hprod_i
    plot(ratetbl[:, i].*2e5, intaltgrid, label=rxnlist[i], color=cullers[c, :], marker=marks[c], ms=7)
    c+=1
end
xscale("log")
legend(bbox_to_anchor=(1.1, 1))
title("Chemical reactions: H production")
xlabel(L"Reaction rate (cm$^{-2}$)")
ylabel("Altitude (km)")
# ylim(25, 75)
savefig(base*"Hproduction.png", bbox_inches="tight")

# D loss and production --------------------------------------------------------
fig = figure(figsize=(8,6))
cullers = get_colors(length(Dloss_i))
c = 1
for i in Dloss_i
    plot(ratetbl[:, i].*2e5, intaltgrid, label=rxnlist[i], color=cullers[c, :], marker=marks[c], ms=7)
    c+=1
end
xscale("log")
legend(bbox_to_anchor=(1.1, 1))
title("Chemical reactions: D loss")
ylim(25, 75)
xlabel(L"Reaction rate (cm$^{-2}$)")
ylabel("Altitude (km)")
savefig(base*"Dloss.png", bbox_inches="tight")

fig = figure(figsize=(8,6))
cullers = get_colors(length(Dprod_i))
c = 1
for i in Dprod_i
    plot(ratetbl[:, i].*2e5, intaltgrid, label=rxnlist[i], color=cullers[c, :], marker=marks[c], ms=7)
    c+=1
end
xscale("log")
legend(bbox_to_anchor=(1.1, 1))
title("Chemical reactions: D production")
#ylim(25, 75)
xlabel(L"Reaction rate (cm$^{-2}$)")
ylabel("Altitude (km)")
savefig(base*"Dproduction.png", bbox_inches="tight")