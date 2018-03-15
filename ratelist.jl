# reactions and multipliers on base rates for deuterium reactions from Yung
# 1988; base rates from this work or Chaffin+ 2017

[[:CO2],[:CO,:O],:JCO2toCOpO],
[[:CO2],[:CO,:O1D],:JCO2toCOpO1D],
[[:O2],[:O,:O],:JO2toOpO],
[[:O2],[:O,:O1D],:JO2toOpO1D],
[[:O3],[:O2,:O],:JO3toO2pO],
[[:O3],[:O2,:O1D],:JO3toO2pO1D],
[[:O3],[:O,:O,:O],:JO3toOpOpO],
[[:H2],[:H,:H],:JH2toHpH],
[[:OH],[:O,:H],:JOHtoOpH],
[[:OH],[:O1D,:H],:JOHtoO1DpH],
[[:HO2],[:OH,:O],:JHO2toOHpH], # other branches should be here, but
                             # have not been measured
# Photolysis of water and heavy water
[[:H2O], [:H,:OH], :JH2OtoHpOH],
[[:H2O], [:H2,:O1D], :JH2OtoH2pO1D],
[[:H2O], [:H,:H,:O], :JH2OtoHpHpO],
[[:HDO], [:H, :OD], :JHDOtoHpOD],
[[:HDO], [:D, :OH], :JHDOtoDpOH],

# Photolysis of hydrogen peroxide and analogues
[[:H2O2], [:OH, :OH], :JH2O2to2OH],
[[:H2O2], [:HO2, :H], :JH2O2toHO2pH],
[[:H2O2], [:H2O, :O1D], :JH2O2toH2OpO1D],
[[:HDO2], [:OH, :OD], :JHDO2toOHpOD],

# recombination of O
[[:O,:O,:M],[:O2,:M],:(1.8*3.0e-33*(300/T)^3.25)],
[[:O,:O2,:N2],[:O3,:N2],:(5e-35*exp(724./T))],
[[:O,:O2,:CO2],[:O3,:CO2],:(2.5*6.0e-34*(300./T)^2.4)],
[[:O,:O3],[:O2,:O2],:(8.0e-12*exp(-2060./T))],
[[:O,:CO,:M],[:CO2,:M],:(2.2e-33*exp(-1780./T))],

# O1D attack
[[:O1D,:O2], [:O,:O2], :(3.2e-11*exp(70./T))],
[[:O1D,:O3], [:O2,:O2], :(1.2e-10)],
[[:O1D,:O3], [:O,:O,:O2], :(1.2e-10)],
[[:O1D, :CO2], [:O,:CO2], :(7.5e-11*exp(115./T))],
## O1D + H2
[[:O1D, :H2],[:H,:OH], :(1.2e-10)],
[[:O1D, :HD], [:H, :OD], :(0.41*1.2e-10)],
[[:O1D, :HD], [:D, :OH], :(0.41*1.2e-10)],
## O1D + H2O
[[:O1D, :H2O],[:OH,:OH], :(1.63e-10*exp(60./T))],
[[:O1D, :HDO], [:OD, :OH], :(1.63e-10*exp(60./T))],

# loss of H2
[[:H2,:O], [:OH,:H], :(6.34e-12*exp(-4000/T))],
### [[:HD,:O], [:OD,:H], :()],
### [[:HD,:O], [:OH,:D], :()],
## OH + H2
[[:OH, :H2], [:H2O, :H], :(7.7E-12*exp(-2100/T))],  # new, KIDA. Old: 9.01e-13*exp(-1526/T)
[[:OD, :H2], [:HDO, :H], :(7.7E-12*exp(-2100/T))], #old: 9.01e-13*exp(-1526/T), same as OH+H2->H2O+H
#[[:OD, :H2], [:H2O, :D], :(0)],
[[:OH, :HD], [:HDO, :H], :((3./20)*7.7E-12*exp(-2100/T))], # old base: 9.01e-13*exp(-1526/T). Thi also has different rate for this rxn
[[:OH, :HD], [:H2O, :D], :((3./20)*7.7E-12*exp(-2100/T))],

# recombination of H
[[:H, :H, :M], [:H2, :M], :(1.5e-29*(T^-1.3))], # NEW - from Yung 88 to match the D version
[[:H, :D, :M], [:HD, :M], :(1.5e-29*(T^-1.3))], # FIX: CHECK THESE RATES??***
### possible deuterated analogues below
[[:H, :H, :CO2], [:H2, :CO2], :(1.6e-32*(298./T)^2.27)],
### [[:H, :D, :CO2], [:HD, :CO2], :()],
[[:H, :OH, :CO2], [:H2O, :CO2], :(1.9*6.8e-31*(300./T)^2)],
### [[:H, :OD, :CO2], [:HDO, :CO2], :()],
### [[:D, :OH, :CO2], [:HDO, :CO2], :()],

## H + HO2
[[:H, :HO2], [:OH, :OH], :(7.2e-11)],
[[:H, :HO2], [:H2, :O2], :(5.6e-12)], # for 245-300K. new rate from KIDA. old: 0.5*6.9e-12. alt new: 1.75e-10*exp(-1030/T) 250-1000K
[[:H, :HO2], [:H2O, :O1D], :(1.6e-12)], # O1D is theoretically mandated
[[:H, :DO2], [:OH, :OD], :(7.2e-11)],
[[:H, :DO2], [:HD, :O2], :(5.6e-12)], # old: 0.5*6.9e-12. same as H+HO2->H2+O2 # FIX THIS!!!!!!!!!!!!!!!!!!!!!!***
[[:H, :DO2], [:HDO, :O1D], :(1.6e-12)],
[[:D, :HO2], [:OH, :OD], :(0.71*7.2e-11)],
[[:D, :HO2], [:HD, :O2], :(0.71*5.6e-12)], # old base: 6.9e-12. Changed to match analogues above
[[:D, :HO2], [:HDO, :O1D], :(0.71*1.6e-12)], # also updated the species to O1D to match H analogue reaction
                                             # change to O if more troubleshooting needed
[[:H, :DO2], [:HO2, :D], :(1e-10/(0.54*exp(890/T)))],
[[:D, :HO2], [:DO2, :H], :(1.0e-10)],

## H + H2O2. possible deuterated analogues below
[[:H, :H2O2], [:HO2,:H2],:(2.8e-12*exp(-1890/T))],
### [[:H, :HDO2], [:DO2,:H2], :()],
### [[:H, :HDO2], [:HO2,:HD], :()],
### [[:H, :HDO2], [:HDO,:OH], :()],
### [[:H, :HDO2], [:H2O,:OD], :()],
[[:H, :H2O2], [:H2O,:OH],:(1.7e-11*exp(-1800/T))],
### [[:D, :H2O2], [:DO2,:H2], :()],
### [[:D, :H2O2], [:HO2,:HD], :()],
### [[:D, :H2O2], [:HDO,:OH], :()],
### [[:D, :H2O2], [:H2O,:OD], :()],

# Interconversion of odd H
## H + O2
[[:H, :O2], [:HO2], threebody(:(2.0*4.4e-32*(T/300)^-1.3),
                          :(7.5e-11*(T/300)^0.2))],
[[:D, :O2], [:DO2], threebody(:(2.0*4.4e-32*(T/300)^-1.3),
                          :(7.5e-11*(T/300)^0.2))],
## H + O3
[[:H, :O3], [:OH, :O2], :(1.4e-10*exp(-470/T))],
[[:D, :O3], [:OD, :O2], :(0.71*1.4e-10*exp(-470/T))],
## O + OH
[[:O, :OH], [:O2, :H], :(2.4e-11*exp(110/T))], # KIDA, for T in 150-500K. old: 1.8e-11*exp(180/T)
[[:O, :OD], [:O2, :D], :(2.4e-11*exp(110/T))], #old: 1.8e-11*exp(180/T). same as O+OH->O2+H
## O + HO2
[[:O, :HO2], [:OH, :O2], :(2.7e-11*exp(224/T))], # KIDA for 220-400K. old: 3.0e-11*exp(200./T)
[[:O, :DO2], [:OD, :O2], :(2.7e-11*exp(224/T))], #old: 3.0e-11*exp(200./T). same as O+HO2 -> OH+O2
## O + H2O2
[[:O, :H2O2], [:OH, :HO2], :(1.4e-12*exp(-2000/T))],
[[:O, :HDO2], [:OD, :HO2], :(0.5*1.4e-12*exp(-2000/T))],
[[:O, :HDO2], [:OH, :DO2], :(0.5*1.4e-12*exp(-2000/T))],
## OH + OH
[[:OH, :OH], [:H2O, :O], :(1.65e-12*(T/300)^1.14*exp(-50/T))], # from UMIST. old: 1.8e-12. KIDA: 6.2e-14*(T/300)^2.62*exp(945/T) for 200-350K
[[:OD, :OH], [:HDO, :O], :(1.65e-12*(T/300)^1.14*exp(-50/T))], #  old: 1.8e-12. same as OH + OH -> H2O + O
[[:OH, :OH], [:H2O2], threebody(:(1.3*6.9e-31*(T/300)^-1.0),:(2.6e-11))],
[[:OD, :OH], [:HDO2], threebody(:(1.3*6.9e-31*(T/300)^-1.0),:(2.6e-11))],
## OH + O3
[[:OH, :O3], [:HO2, :O2], :(1.7e-12*exp(-940/T))],
[[:OD, :O3], [:DO2, :O2], :(1.7e-12*exp(-940/T))],
## OH + HO2
[[:OH, :HO2], [:H2O, :O2], :(4.8e-11*exp(250/T))],
[[:OH, :DO2], [:HDO, :O2], :(4.8e-11*exp(250/T))],
[[:OD, :HO2], [:HDO, :O2], :(4.8e-11*exp(250/T))],
## OH + H2O2
[[:OH, :H2O2], [:H2O, :HO2], :(5.26e-12*exp(-307/T))], #from UMIST. old: 1.8e-12. KIDA: 2.9e-12*exp(-160/T) for 240-460K
[[:OD, :H2O2], [:HDO, :HO2], :(5.26e-12*exp(-307/T))], # old: 1.8e-12.
[[:OH, :HDO2], [:HDO, :HO2], :(0.5*5.26e-12*exp(-307/T))], # old: 0.5*1.8e-12
[[:OH, :HDO2], [:H2O, :DO2], :(0.5*5.26e-12*exp(-307/T))], # see previous.
#[[:OD, :H2O2], [:H2O, :DO2], :(0)],
## HO2 + O3
[[:HO2, :O3], [:OH, :O2, :O2], :(2.03e-16*(T/300)^4.57*exp(693/T))], # from KIDA for 250-340K. old: 1.0e-14*exp(-490/T)
[[:DO2, :O3], [:OD, :O2, :O2], :(2.03e-16*(T/300)^4.57*exp(693/T))],  # old: 1.0e-14*exp(-490/T)
## HO2 + HO2
[[:HO2, :HO2], [:H2O2, :O2], :(2.2e-13*exp(600/T))], # from KIDA for 230-420K. old: 3.0e-13*exp(460/T). UMIST: 1e-12
[[:DO2, :HO2], [:HDO2, :O2], :(2.2e-13*exp(600/T))], # old: 3.0e-13*exp(460/T)
## HO2 + HO2 + M
[[:HO2, :HO2, :M], [:H2O2, :O2, :M],:(2*2.1e-33*exp(920/T))],
### [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :()],

## OH + D or OD + H (no non-deuterated analogues)
[[:OD, :H], [:OH, :D], :(3.3e-9*(T^-0.63)/(0.72*exp(717/T)))], # Thi: 1.26e-10*(T/300)^-0.63*exp(-717/T)
[[:OH, :D], [:OD, :H], :(3.3e-9*T^-0.63)], # Thi: 9.07e-11*(T/300)^-0.63

# CO2 recombination due to odd H (with HOCO intermediate)
## straight to CO2
[[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13*(T/300)^0.6),:(2.1e9*(T/300)^6.1))], # from Sander
[[:OD, :CO], [:CO2, :D], threebodyca(:(1.5e-13*(T/300)^0.6),:(2.1e9*(T/300)^6.1))],
### possible deuterated analogues below
[[:OH,:CO], [:HOCO], threebody(:(5.9e-33*(T/300)^-1.4),:(1.1e-12*(T/300)^1.3))],
### [[:OD,:CO], [:DOCO], threebody(:(),:()], # WOULD ADD NEW SPECIES
[[:HOCO, :O2], [:HO2,:CO2], :(2.0e-12)], # not found in databases
### [[:DOCO, :O2], [:DO2,:CO2], :()], # WOULD ADD NEW SPECIES

# CO2+ attack on molecular hydrogen
[[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # UMIST: results in HCO2pl + H with rate 9.5e-10
[[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5)*8.7e-10)],
