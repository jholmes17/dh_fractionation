################################################################################
# converge_new_file.jl
# TYPE: MAIN (main script)
# WHICH: Equilibrium and perturbation experiments
# DESCRIPTION: does initial convergence for a photochemistry experiment of the
# Martian atmosphere.
# 
# Eryn Cangi
# 2018-2019
# this script should be used with converged_standarwater.h5 and called before
# run_and_plot.jl. 
# Currently tested for Julia: 0.7
################################################################################

###ALL distance/area/volume units must be in CM!!!!

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra

################################################################################
############################# FUNCTION DEFINITIONS #############################
################################################################################

# misc. utility functions ======================================================
function fluxsymbol(x)
    Symbol(string("f",string(x)))
end

function get_ncurrent(readfile)
    #=
    Retrieves the matrix of species concentrations by altitude from an HDF5
    file containing a converged atmosphere.
    =#
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()
    
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    return n_current
end

function write_ncurrent(n_current, filename)
    #=
    Writes out the current state of species concentrations by altitude to a file
    (for converged atmosphere). 
    =# 
    n_current_mat = Array{Float64}(undef, length(alt)-2, 
                                   length(collect(keys(n_current))));
    for ispecies in [1:length(collect(keys(n_current)));]
        for ialt in [1:length(alt)-2;]
            n_current_mat[ialt, ispecies] = n_current[collect(keys(n_current))[ispecies]][ialt]
        end
    end
    h5write(filename,"n_current/n_current_mat",n_current_mat)
    h5write(filename,"n_current/alt",alt)
    h5write(filename,"n_current/species",map(string, collect(keys(n_current))))
end

function readandskip(f, delim::Char, T::Type; skipstart=0)
    #= 
    function to read in data from a file, skipping zero or more lines at the 
    beginning.
    =# 
    f = open(f,"r")
    if skipstart>0
        for i in [1:skipstart;]
            readline(f)
        end
    end
    f = readdlm(f, delim, T)
end

# Density functions ============================================================

function n_tot(n_current, z)
    #= get the total number density at a given altitude =#
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

function meanmass(n_current, z)
    #= find the mean molecular mass at a given altitude =#
    thisaltindex = n_alt_index[z]
    c = [n_current[sp][thisaltindex] for sp in specieslist]
    m = [speciesmolmasslist[sp] for sp in specieslist]
    return sum(c.*m)/sum(c)
end

# transport/scale height =======================================================

function scaleH(z, T::Float64, mm::Real)
    #= Computes the scale height of the atmosphere for the mean molar mass =#
    return boltzmannK*T/(mm*mH*marsM*bigG)*(((z+radiusM)*1e-2)^2)*1e2
    # constants are in MKS. Convert to m and back to cm.
end

function scaleH(z, species::Symbol)
    #=
    computes scale height of the atmosphere for a specific species
    =#
    T=Temp(z)
    mm = speciesmolmasslist[species]
    scaleH(z, T, mm)
end

function scaleH(z, T::Float64, species::Symbol)
    #=
    computes scale height of the atmosphere for a specific species
    =#
    mm = speciesmolmasslist[species]
    scaleH(z, T, mm)
end

function scaleH(z, T::Float64, n_current)
    #= Computes the scale height of the atmosphere for the mean molar mass =#
    mm = meanmass(n_current, z)
    scaleH(z, T, mm)
end

# transport equations ==========================================================

#=
    at each level of the atmosphere, density can be transferred up or
    down with a specified rate coefficient.

                             | n_i+1
          ^  tspecies_i_up   v tspecies_i+1_down
      n_i |
          v tspecies_i_down  ^ tspecies_i-1_up
                             | n_i-1

    the flux at each cell boundary is the sum of the upward flux from
    the cell below and the downward flux of the cell above. These fluxes
    are determined using flux coefficients that come from the diffusion
    equation. Care must be taken at the upper and lower boundary so that
    tspecies_top_up and tspecies_bottom_down properly reflect the
    boundary conditions of the atmosphere.

    This is handled in the code with the population of the appropriate
    reactions, with variable rate coeffecients that are populated
    between each timestep (similar to the way photolysis rates are
    included). We need to make reactions at each interior altitude
    level:
             n_i -> n_i+1  tspecies_i_up
             n_i -> n_i-1  tspecies_i_down

    At the upper and lower boundary we omit the species on the RHS, so
    that these reactions are potentially non-conservative:

            n_top    -> NULL  tspecies_top_up
            n_bottom -> NULL  tspecies_bottom_down

    These coefficients then describe the diffusion velocity at the top
    and bottom of the atmosphere.
=#

function fluxcoefs(z, dz, Kv, Dv, Tv, Hsv, H0v, species) #ntv,
    #= 
    function to generate coefficients of the transport network.
    z: altitude
    dz: altitude layer thickness ("resolution")
    Kv: eddy diffusion coefficient
    Dv: molecular diffusion coefficient
    Tv: temperature
    Hsv: scale height by species
    H0v: mean atmospheric scale height
    species: species symbol 

    v refers to "vector"
    u refers to "upper" (the upper parcel)
    l refers to "lower" (the lower parcel)
    =#
    Dl = (Dv[1] + Dv[2])/2.0
    Kl = (Kv[1] + Kv[2])/2.0
    Tl = (Tv[1] + Tv[2])/2.0
    dTdzl = (Tv[2] - Tv[1])/dz
    Hsl = (Hsv[1] + Hsv[2])/2.0
    H0l = (H0v[1] + H0v[2])/2.0

    # we have two flux terms to combine:
    sumeddyl = (Dl+Kl)/dz/dz
    gravthermall = (Dl*(1/Hsl + (1+thermaldiff(species))/Tl*dTdzl) +
                    Kl*(1/H0l + 1/Tl*dTdzl))/(2*dz)

    Du = (Dv[2] + Dv[3])/2.0
    Ku = (Kv[2] + Kv[3])/2.0
    Tu = (Tv[2] + Tv[3])/2.0
    dTdzu = (Tv[3] - Tv[2])/dz
    Hsu = (Hsv[2] + Hsv[3])/2.0
    H0u = (H0v[2] + H0v[3])/2.0

    # we have two flux terms to combine:
    sumeddyu = (Du+Ku)/dz/dz
    gravthermalu = (Du*(1/Hsu + (1 + thermaldiff(species))/Tu*dTdzu) +
                  Ku*(1/H0u + 1/Tu*dTdzu))/(2*dz)

    # this results in the following coupling coefficients:
    return [sumeddyl+gravthermall, # down
            sumeddyu-gravthermalu] # up
end

function fluxcoefs(z, dz, species, n_current)
    # overload to generate the coefficients if they are not supplied
    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)
    Tp = Temp(z+dz)
    T0 = Temp(z)
    Tm = Temp(z-dz)
    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = scaleH(z+dz, species)
    Hs0 = scaleH(z, species)
    Hsm = scaleH(z-dz, species)
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)
end

function lower_up(z, dz, species, n_current)
    #= 
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one above. Variables ending in p refer to
    the layer above ("plus"), in 0 refer to the current layer at altitude z,
    and ending in m refer to the layer below ("minus") 

    z: altitude
    dz: altitude layer thickness ("resolution")
    species: species symbol
    n_current: current atmospheric state in a matrix by altitude and species.
    =#
    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = 1
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z,nt0)
    Km = 1
    Tp = Temp(z+dz)
    T0 = Temp(z)
    Tm = 1
    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = 1
    Hsp = scaleH(z+dz, species)
    Hs0 = scaleH(z,species)
    Hsm = 1
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z,T0, n_current)
    H0m = 1

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[2]
end

function upper_down(z, dz, species, n_current)
    #= 
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one above. Variables ending in p refer to
    the layer above ("plus"), in 0 refer to the current layer at altitude z,
    and ending in m refer to the layer below ("minus") 

    z: altitude
    dz: altitude layer thickness ("resolution")
    species: species symbol
    n_current: current atmospheric state in a matrix by altitude and species.
    =#
    ntp = 1
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = 1
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)
    Tp = 1
    T0 = Temp(z)
    Tm = Temp(z-dz)
    Dp = 1
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = 1
    Hs0 = scaleH(z, species)
    Hsm = scaleH(z-dz, species)
    H0p = 1
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[1]
end

function boundaryconditions(species, dz, n_current)
    #= returns the symbolic transport coefficients that encode the
    boundary conditions for the null-pointing equations

    n_1->NULL t_lower_bc
    n_f->NULL t_upper_bc

    this defines two additional symbols for each species that need
    to be resolved in the function call macro:
                tspecies_lower_up and
                tspecies_upper_down
    these are found by passing the appropriate values to fluxcoefs
    and selecting the correct output
    =#

    bcs = speciesbcs(species)
    if issubset([species],notransportspecies)
        bcs = ["f" 0.; "f" 0.]
    end

    #= first element returned corresponds to lower BC, second to upper
    BC transport rate. Within each element, the two rates correspond
    to the two equations
    n_b  -> NULL (first rate, depends on species concentration)
    NULL -> n_b  (second rate, independent of species concentration =#
    bcvec = Float64[0 0;0 0]

    # LOWER
    if bcs[1, 1] == "n"
        bcvec[1,:]=[fluxcoefs(alt[2], dz, species, n_current)[1],
                    lower_up(alt[1], dz, species, n_current)*bcs[1, 2]]
    elseif bcs[1, 1] == "f"
        bcvec[1,:] = [0.0, bcs[1, 2]/dz]
    elseif bcs[1, 1] == "v"
        bcvec[1,:] = [bcs[1, 2]/dz, 0.0]
    else
        throw("Improper lower boundary condition!")
    end

    # UPPER
    if bcs[2, 1] == "n"
        bcvec[2,:] = [fluxcoefs(alt[end-1],dz, species, n_current)[2],
                    upper_down(alt[end],dz, species, n_current)*bcs[1, 2]]
    elseif bcs[2, 1] == "f"
            bcvec[2,:] = [0.0,-bcs[2, 2]/dz]
    elseif bcs[2, 1] == "v"
        bcvec[2,:] = [bcs[2, 2]/dz, 0.0]
    else
        throw("Improper upper boundary condition!")
    end

    # return the bc vec
    return bcvec
end

# chemistry equations ==========================================================

# note: this whole section is basically witchcraft that Mike wrote and I've never
# really looked into it in detail. I just trust it works. --Eryn

function loss_equations(network, species)
    #=  
    given a network of equations in the form of reactionnet above, this
    function returns the LHS (reactants) and rate coefficient for all
    reactions where the supplied species is consumed, in the form of an array
    where each entry is of the form [reactants, rate] 
    =#

    # get list of all chemical reactions species participates in:
    speciespos = getpos(network, species)
    # find pos where species is on LHS but not RHS:
    lhspos = map(x->x[1],  # we only need the reaction number
               map(x->speciespos[x],  # select the appropriate reactions
                   findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1],  # we only need the reaction number
               map(x->speciespos[x],  # select the appropriate reactions
                   findall(x->x[2]==2, speciespos)))
    for i in intersect(lhspos, rhspos)
        lhspos = deletefirst(lhspos, i)
    end

    # get the products and rate coefficient for the identified reactions.
    losseqns=map(x->vcat(Any[network[x][1]...,network[x][3]]), lhspos)
    # automatically finds a species where it occurs twice on the LHS
end

function loss_rate(network, species)
    #= return a symbolic expression for the loss rate of species in the
    supplied reaction network. Format is a symbolic expression containing a sum
    of reactants * rate. =#
    leqn=loss_equations(network, species) # get the equations
    lval=:(+($( # and add the products together
               map(x->:(*($(x...))) # take the product of the
                                    # concentrations and coefficients
                                    # for each reaction
                   ,leqn)...)))
end

function production_equations(network, species)
    #= given a network of equations in the form of reactionnet above, this
    function returns the LHS (reactants) and rate coefficient for all
    reactions where the supplied species is a product, in the form of an array
    where each entry is of the form [reactants, rate] =#

    speciespos = getpos(network, species)#list of all reactions where species is produced
    # find pos where species is on RHS but not LHS
    lhspos = map(x->x[1], # we only need the reaction number
               map(x->speciespos[x], #select the appropriate reactions
                   findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1], # we only need the reaction number
               map(x->speciespos[x], # select the appropriate reactions
                   findall(x->x[2]==2, speciespos)))
    for i in intersect(rhspos, lhspos)
        rhspos = deletefirst(rhspos, i)
    end

    # get the products and rate coefficient for the identified reactions.
    prodeqns = map(x->vcat(Any[network[x][1]...,network[x][3]]),
                 # automatically finds and counts duplicate
                 # production for each molecule produced
                 rhspos)

    return prodeqns
end

function production_rate(network, species)
    #= 
    return a symbolic expression for the loss rate of species in the
    supplied reaction network.
    =#

    # get the reactants and rate coefficients
    peqn = production_equations(network, species)

    # add up and take the product of each set of reactants and coeffecient
    pval = :(+ ( $(map(x -> :(*($(x...))), peqn) ...) ))
end

function chemical_jacobian(chemnetwork, transportnetwork, specieslist, dspecieslist)
    #= 
    Compute the symbolic chemical jacobian of a supplied reaction network
    for the specified list of species. Returns three arrays suitable for
    constructing a sparse matrix: lists of the first and second indices
    and the symbolic value to place at that index.
    =#

    # set up output vectors: indices and values
    ivec = Int64[] # list of first indices (corresponding to the species being produced and lost)
    jvec = Int64[] # list of second indices (corresponding to the derivative being taken)
    tvec = Any[] # list of the symbolic values corresponding to the jacobian

    nspecies = length(specieslist)
    ndspecies = length(dspecieslist)
    for i in 1:nspecies #for each species
        ispecies = specieslist[i]
        # get the production and loss equations
        peqn = []
        leqn = []
        if issubset([ispecies],chemspecies)
            peqn = [peqn; production_equations(chemnetwork, ispecies)] #***
            leqn = [leqn; loss_equations(chemnetwork, ispecies)]
        end
        if issubset([ispecies],transportspecies)
            peqn = [peqn; production_equations(transportnetwork, ispecies)]
            leqn = [leqn; loss_equations(transportnetwork, ispecies)]
        end
        for j in 1:ndspecies #now take the derivative with resp`ect to the other species
            jspecies = dspecieslist[j]
            #= find the places where the production rates depend on
            jspecies, and return the list rates with the first
            occurrance of jspecies deleted. (Note: this seamlessly
            deals with multiple copies of a species on either side of
            an equation, because it is found twice wherever it lives) =#
            ppos = map(x->deletefirst(peqn[x[1]],jspecies),getpos(peqn, jspecies))
            lpos = map(x->deletefirst(leqn[x[1]],jspecies),getpos(leqn, jspecies))
            if length(ppos)+length(lpos)>0 #if there is a dependence
                #make note of where this dependency exists
                append!(ivec,[i])
                append!(jvec,[j])
                #= smash the production and loss rates together,
                multiplying for each distinct equation, adding
                together the production and loss seperately, and
                subtracting loss from production. =#
                if length(ppos)==0
                    lval = :(+($(map(x->:(*($(x...))),lpos)...)))
                    tval = :(-($lval))
                elseif length(lpos)==0
                    pval = :(+($(map(x->:(*($(x...))),ppos)...)))
                    tval = :(+($pval))
                else
                    pval = :(+($(map(x->:(*($(x...))),ppos)...)))
                    lval = :(+($(map(x->:(*($(x...))),lpos)...)))
                    tval = :(-($pval,$lval))
                end
                # attach the symbolic expression to the return values
                append!(tvec,[tval])
            end
        end
    end
    return (ivec, jvec, Expr(:vcat, tvec...))
end

function getrate(chemnet, transportnet, species)
    #=
    Creates a symbolic expression for the rate at which a given species is
    either produced or lost. Production is from chemical reaction yields or
    entry from other atmospheric layers. Loss is due to consumption in reactions
     or migration to other layers.
    =#
    rate = :(0.0)
    if issubset([species],chemspecies)
        rate = :($rate
               + $(production_rate(chemnet, species))
               - $(      loss_rate(chemnet, species)))
    end
    if issubset([species],transportspecies)
        rate = :($rate
               + $(production_rate(transportnet, species))
               - $(      loss_rate(transportnet, species)))
    end

    return rate
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
    return theserates
end

function getflux(n_current, dz, species)
    #=
    Returns a 1D array of fluxes in and out of a given altitude level for a given species
    =#
    thesecoefs = [fluxcoefs(a, dz, species, n_current) for a in alt[2:end-1]]
    thesebcs = boundaryconditions(species, dz, n_current)

    thesefluxes = fill(convert(Float64, NaN),length(intaltgrid))

    thesefluxes[1] = (-(n_current[species][2]*thesecoefs[2][1]
                      -n_current[species][1]*thesecoefs[1][2])
                    +(-n_current[species][1]*thesebcs[1, 1]
                      +thesebcs[1, 2]))/2.0
    for ialt in 2:length(intaltgrid)-1
        thesefluxes[ialt] = (-(n_current[species][ialt+1]*thesecoefs[ialt+1][1]
                             - n_current[species][ialt]*thesecoefs[ialt][2])
                             + (-n_current[species][ialt]*thesecoefs[ialt][1]
                             + n_current[species][ialt-1]*thesecoefs[ialt-1][2]))/2.0
    end
    thesefluxes[end] = (-(thesebcs[2, 2]
                        - n_current[species][end]*thesebcs[2, 1])
                        + (-n_current[species][end]*thesecoefs[end][1]
                        + n_current[species][end-1]*thesecoefs[end-1][2]))/2.0
    return dz*thesefluxes
end

function fluxes(n_current, dz)
    thesefluxes = fill(convert(Float64, NaN),(length(intaltgrid),length(specieslist)))
    for isp in 1:length(specieslist)
        thesefluxes[:,isp] = getflux(n_current, dz, specieslist[isp])
    end
    thesefluxes
end

function ratefn(nthis, inactive, Jrates, T, M, tup, tdown, tlower, tupper)
    # at each altitude, get the appropriate group of concentrations,
    # coefficients, and rates to pass to ratefn_local
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive,(length(inactivespecies),length(intaltgrid)))

    returnrates = zero(nthismat)

    # fill the first altitude entry with information for all species
    returnrates[:,1] = ratefn_local([nthismat[:,1]; nthismat[:,2];
                                    fill(1.0, length(activespecies));
                                    inactivemat[:,1]; Jrates[:,1]; T[1]; M[1];
                                    tup[:,1]; tlower[:,1]; tdown[:,2];
                                    tlower[:,2]]...)

    # iterate through other altitudes except the last level, filling the info in
    for ialt in 2:(length(intaltgrid)-1)
        returnrates[:,ialt] = ratefn_local([nthismat[:,ialt];
                                          nthismat[:,ialt+1];
                                          nthismat[:,ialt-1];
                                          inactivemat[:,ialt];
                                          Jrates[:,ialt];
                                          T[ialt]; M[ialt];
                                          tup[:,ialt]; tdown[:,ialt];
                                          tdown[:,ialt+1]; tup[:,ialt-1]]...)
    end

    # fill in the last level of altitude (200 km)
    returnrates[:,end] = ratefn_local([nthismat[:,end];
                                       fill(1.0, length(activespecies));
                                       nthismat[:,end-1];
                                       inactivemat[:,end];
                                       Jrates[:,end];
                                       T[end]; M[end];
                                       tupper[:,1]; tdown[:,end];
                                       tupper[:,2]; tup[:,end-1]]...)
    return [returnrates...;]
end

function chemJmat(nthis, inactive, Jrates, T, M, tup, tdown, tlower, tupper, dt)
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive, (length(inactivespecies), length(intaltgrid)))
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]

    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,1]; nthismat[:,2];
                                                  fill(1.0, length(activespecies));
                                                  inactivemat[:,1]; Jrates[:,1];
                                                  T[1]; M[1]; tup[:,1]; tlower[:,1];
                                                  tdown[:,2]; tlower[:,2];dt]...)
    #add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])
    #and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(activespecies))
    append!(chemJval, tcupper[3])

    for ialt in 2:(length(intaltgrid)-1)
        (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,ialt];
                                                      nthismat[:,ialt+1];
                                                      nthismat[:,ialt-1];
                                                      inactivemat[:,ialt];
                                                      Jrates[:,ialt]; T[ialt]; 
                                                      M[ialt]; tup[:,ialt]; 
                                                      tdown[:,ialt];
                                                      tdown[:,ialt+1];
                                                      tup[:,ialt-1]; dt]...)
        # add the influence of the local densities
        append!(chemJi, tclocal[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclocal[2].+(ialt-1)*length(activespecies))
        append!(chemJval, tclocal[3])
        # and the upper densities
        append!(chemJi, tcupper[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tcupper[2].+(ialt  )*length(activespecies))
        append!(chemJval, tcupper[3])
        # and the lower densities
        append!(chemJi, tclower[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclower[2].+(ialt-2)*length(activespecies))
        append!(chemJval, tclower[3])
    end

    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,end];
                                              fill(1.0, length(activespecies));
                                              nthismat[:,end-1];
                                              inactivemat[:,end];
                                              Jrates[:,end];
                                              T[end]; M[end];
                                              tupper[:,1]; tdown[:,end];
                                              tupper[:,2]; tup[:,end-1];dt]...)
    # add the influence of the local densities
    append!(chemJi, tclocal[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclocal[2].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJval, tclocal[3])
    # and the lower densities
    append!(chemJi, tclower[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclower[2].+(length(intaltgrid)-2)*length(activespecies))
    append!(chemJval, tclower[3])

    # make sure to add 1's along the diagonal
    append!(chemJi,[1:length(nthis);])
    append!(chemJj,[1:length(nthis);])
    append!(chemJval, fill(1.0, length(nthis)))

    sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +);
end

# auxiliary functions ==========================================================

function getpos(array, test::Function, n=Any[])
    #= 
    this function searches through an arbitrarily structured array
    finding elements that match the test function supplied, and returns a
    one-dimensional array of the indicies of these elements. 
    =#
    if !isa(array, Array)
        test(array) ? Any[n] : Any[]
    else
        vcat([ getpos(array[i], test, Any[n...,i]) for i=1:size(array)[1] ]...)
    end
end

function getpos(array, value)
    #= 
    overloading getpos for the most common use case, finding indicies
    corresponding to elements in array that match value. 
    =#
    getpos(array, x->x==value)
end

function deletefirst(A, v)
    # Returns list A with its first element equal to v removed.
    # index = findfirst(A, v) # valid in Julia 0.6, but the order switches in 0.7 ***
    index = something(findfirst(isequal(v), A), 0)
    keep = setdiff([1:length(A);],index)
    return A[keep]
end

# plotting functions ===========================================================
function plotatm(n_current)
    clf()
    # selectspecies = [:O3, :H2O2, :HDO2, :DO2, :HO2, :HD, :H2, :H, :D, :H2O, 
    #                  :HDO, :O, :O2, :CO2]
    for sp in fullspecieslist
        plot(n_current[sp], alt[2:end-1]/1e5, color = speciescolor[sp],
             linewidth=2, label=sp, linestyle=speciesstyle[sp], zorder=1)
    end
    ylim(0, 200)
    xscale("log")
    xlim(1e-15, 1e18)
    xlabel(L"Species concentration (cm$^{-3}$)")
    ylabel("Altitude [km]")
    legend(bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
end

function plotatm()
    plotatm(n_current)
end

# diffusion functions ==========================================================
function Keddy(n_current, z)
    #=
    eddy diffusion coefficient, stolen from Krasnopolsky (1993).
    Scales as the inverse sqrt of atmospheric number density

    n_current: dictionary representing current simulation state.
    z: some altitude in cm.
    =#
    z <= 60.e5 ? 10^6 : 2e13/sqrt(n_tot(n_current, z))
end

function Keddy(z::Real, nt::Real)
    #=
    eddy diffusion coefficient, stolen from Krasnopolsky (1993).
    Scales as the inverse sqrt of atmospheric number density

    z: some altitude in cm.
    nt: number total of species at this altitude (I think)
    =#
    z <= 60.e5 ? 10^6 : 2e13/sqrt(nt)
end

function Dcoef(T, n::Real, species::Symbol)
    #=
    Calculates molecular diffusion coefficient for a particular slice of the
    atmosphere using D = AT^s/n, from Banks and Kockarts Aeronomy, part B,
    pg 41, eqn 15.30 and table 15.2 footnote

    T: temperature
    n: number of molecules (all species) this altitude
    species: whichever species we are calculating for
    =#
    dparms = diffparams(species)
    return dparms[1]*1e17*T^(dparms[2])/n
end

# boundary condition functions =================================================

function effusion_velocity(Texo::Float64, m::Float64)
    #=
    Returns effusion velocity for a species.
    Texo: temperature of the exobase (upper boundary) in K
    m: mass of one molecule of species in amu
    =#
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+zmax))
    vth = sqrt(2*boltzmannK*Texo/(m*mH))
    v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)
    return v
end

function speciesbcs(species)
    get(speciesbclist,
        species,
        ["f" 0.; "f" 0.])
end

# main routine functions =======================================================
function update_Jrates!(n_current::Dict{Symbol, Array{Float64, 1}})
    #=
    this function updates the photolysis rates stored in n_current to
    reflect the altitude distribution of absorbing species
    =#
    #    global solarabs::Array{Array{Float64, 1},1}

    # Initialize an array, length=num of altitude levels - 2.
    # Each sub-array is an array of length 2000, corresponding to 2000 wavelengths.
    solarabs = Array{Array{Float64}}(undef, length(alt)-2)
    for i in range(1, length=length(alt)-2)
        solarabs[i] = zeros(Float64, 2000)
    end

    nalt = size(solarabs, 1)
    nlambda = size(solarabs[1],1)

    for jspecies in Jratelist
        species = absorber[jspecies]

        jcolumn = 0.
        for ialt in [nalt:-1:1;]
            #get the vertical column of the absorbing constituient
            jcolumn += n_current[species][ialt]*dz
            # if jspecies==:JO2toOpO
            #     println(string("At alt = ",alt[ialt+1],
            #                    ", n_",species," = ",n_current[species][ialt],
            #                    ", jcolumn = ",jcolumn))
            #     println("and solarabs[ialt] is $(solarabs[ialt]) before we do axpy")
            #     readline(STDIN)
            # end

            # add the total extinction to solarabs:
            # multiplies air column density at all wavelengths by crosssection
            # to get optical depth
            BLAS.axpy!(nlambda, jcolumn, crosssection[jspecies][ialt+1], 1,
                       solarabs[ialt],1)
        end
    end

    #solarabs now records the total optical depth of the atmosphere at
    #each wavelength and altitude

    # actinic flux at each wavelength is solar flux diminished by total
    # optical depth
    for ialt in [1:nalt;]
        solarabs[ialt] = solarflux[:,2].*exp.(-solarabs[ialt])
    end

    #each species absorbs according to its cross section at each
    #altitude times the actinic flux.
    for j in Jratelist
        for ialt in [1:nalt;]
            n_current[j][ialt] = BLAS.dot(nlambda, solarabs[ialt], 1,
                                          crosssection[j][ialt+1], 1)
        end
       # this section for testing sensitivity to J rates
        # if contains(string(j), "H2O") | contains(string(j), "HDO")
        # if contains(string(j), "CO2toCOpO")
           # n_current[j] = n_current[j] ./ 10
        # end
    end
end

function timeupdate(mytime)
    for i = 1:15
        plotatm()
        update!(n_current, mytime)
    end
end

function next_timestep(nstart::Array{Float64, 1}, nthis::Array{Float64, 1},
                       inactive::Array{Float64, 1}, Jrates::Array{Float64, 2},
                       T::Array{Float64, 1}, M::Array{Float64, 1},
                       tup::Array{Float64, 2}, tdown::Array{Float64, 2},
                       tlower::Array{Float64, 2}, tupper::Array{Float64, 2},
                       dt::Float64)
    #=
        moves to the next timestep using Newton's method on the linearized
        coupled transport and chemical reaction network.
    =#
    eps = 1.0 # ensure at least one iteration
    iter = 0
    while eps>1e-8
        nold = deepcopy(nthis)

        # stuff concentrations into update function and jacobian
        fval = nthis - nstart - dt*ratefn(nthis, inactive, Jrates, T, M, tup,
                                          tdown, tlower, tupper)
        updatemat = chemJmat(nthis, inactive, Jrates, T, M, tup, tdown, tlower,
                             tupper, dt)

        # update
        nthis = nthis - (updatemat \ fval)
        # check relative size of update
        eps = maximum(abs.(nthis-nold)./nold)
        iter += 1
        if iter>1e3; throw("too many iterations in next_timestep!"); end;
    end
    return nthis
end

@everywhere function update!(n_current::Dict{Symbol, Array{Float64, 1}},dt)
    # update n_current using the coupled reaction network, moving to
    # the next timestep

    #set auxiliary (not solved for in chemistry) species values, photolysis rates
    inactive = deepcopy(Float64[[n_current[sp][ialt] for sp in inactivespecies, ialt in 1:length(intaltgrid)]...])
    Jrates = deepcopy(Float64[n_current[sp][ialt] for sp in Jratelist, ialt in 1:length(intaltgrid)])

    # extract concentrations
    nstart = deepcopy([[n_current[sp][ialt] for sp in activespecies, ialt in 1:length(intaltgrid)]...])
    M = sum([n_current[sp] for sp in fullspecieslist])

    # set temperatures
    T = Float64[Temp(a) for a in alt[2:end-1]]  

    # take initial guess
    nthis = deepcopy(nstart)

    # get the transport rates
    tup = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current)[2] for sp in specieslist, a in alt[2:end-1]]
    tdown = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current)[1] for sp in specieslist, a in alt[2:end-1]]

    # put the lower bcs and upper bcs in separate arrays; but they are not the
    # right shape! TODO: make this more efficient?
    tlower_temp = [boundaryconditions(sp, dz, n_current)[1,:] for sp in specieslist]
    tupper_temp = [boundaryconditions(sp, dz, n_current)[2,:] for sp in specieslist]

    # reshape tlower and tupper into 2x2 arrays
    tlower = zeros(Float64, length(tlower_temp), 2)
    tupper = zeros(Float64, length(tupper_temp), 2)

    # tlower_temp & tupper_temp have same length; OK to use lower for the range
    for r in range(1, length=length(tlower_temp))
        tlower[r, :] = tlower_temp[r]
        tupper[r, :] = tupper_temp[r]
    end

    # update to next timestep
    nthis = next_timestep(nstart, nthis, inactive, Jrates, T, M, tup, tdown,
                          tlower, tupper, dt)
    nthismat = reshape(nthis,(length(activespecies),length(intaltgrid)))

    # write found values out to n_current
    for s in 1:length(activespecies)
        for ia in 1:length(intaltgrid)
            tn = nthismat[s, ia]
            n_current[activespecies[s]][ia] = tn > 0. ? tn : 0.
        end
    end

    update_Jrates!(n_current)
end


################################################################################
################################## MAIN SETUP ##################################
################################################################################

# Set up simulation files and experiment type ==================================
# Location of main scripts and convergence files. Use 2nd line if on laptop
# scriptdir = "/data/GDrive-CU/Research/chaffincode-working/"
scriptdir = "/home/emc/GDrive-CU/Research/chaffincode-working/"

# Storage location for experiment results
# experimentdir = "/data/GDrive-CU/Research/Results/"
experimentdir = "/home/emc/GDrive-CU/Research/Results/"#/VarWaterTemp/"

# get command line arguments for experiment type. format:
# temp Tsurf Ttropo Texo; water <mixing ratio>; dh <multiplier>; Oflux <cm^-2s^-1>
# examples: temp 190 110 200; water 1e-3; dh 8; Oflux 1.2e8
args = Any[ARGS[i] for i in 1:1:length(ARGS)]

# enter the arguments as extension to the filename (FNext); 
if args[1]=="temp"
    FNext = "temp_$(args[2])_$(args[3])_$(args[4])"
elseif args[1]=="water"
    FNext = "water_$(args[2])"
elseif args[1]=="dh"
    FNext = "dh_$(args[2])"
elseif args[1]=="Oflux"
    FNext = "Oflux_$(args[2])"
end

# Set up the folder if it doesn't exist
println("ALERT: running sim for $(FNext)")
println("Checking for existence of $(FNext) folder in results")
dircontents = readdir(experimentdir)
if FNext in dircontents
    println("Folder $(FNext) exists")
else
    mkdir(experimentdir*FNext)
    println("Created folder ", FNext)
end

# Set up the converged file to read from and load the simulation state at init.
readfile = scriptdir*"converged_standardwater.h5"
println("ALERT: Using file: ", readfile)
const alt = h5read(readfile,"n_current/alt")
n_current = get_ncurrent(readfile)
n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])

# set SVP to be fixed or variable with temperature
if args[1] == "temp"
    fix_SVP = true
    println("ALERT: SVP will be fixed to that of the mean temperature profile ")
else
    fix_SVP = false  # temp doesn't vary in other experiments, so false for simplicity
end

# converts the strings that are numeric into numbers, needed to use in functions
for i in 2:1:length(args)
    args[i] = parse(Float64, args[i])
end

# Plot styles ==================================================================
rcParams = PyCall.PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 22
rcParams["axes.labelsize"]= 24
rcParams["xtick.labelsize"] = 22
rcParams["ytick.labelsize"] = 22

 # Style dictionaries for plotting
 speciescolor = Dict( # H group
                      :H => "#ff0000", :D => "#ff0000", # red
                      :H2 => "#e526d7", :HD =>  "#e526d7", # dark pink/magenta

                      # hydroxides
                      :OH => "#7700d5", :OD => "#7700d5", # purple

                      # water group (roughly, I ain't a chemist)
                      :H2O => "#0083dc", :HDO => "#0083dc", # cornflower blue
                      :H2O2 => "#0000ff", :HDO2 => "#0000ff", # true blue
                      :HO2 => "#046868", :DO2 => "#046868",  # dark teal

                      # O group
                      :O1D => "#808000", # olive
                      :O => "#1a6115",   # forest green
                      :O2 => "#15da09",  # kelly/grass green
                      :O3 => "#269e56",  # light green

                      # CO group
                      :CO2 => "#d18564",   # dark peach
                      :CO2pl => "#614215", # brown
                      :CO => "#ff6600",    # orange
                      :HOCO => "#e8ba8c", :DOCO => "#e8ba8c",  #tannish

                      # nonreactants
                      :Ar => "#808080", :N2 => "#cccccc",);  # grays

 speciesstyle = Dict( # H group
                      :H => "-", :D => "--",
                      :H2 => "-",  :HD => "--",
                      # hydroxides
                      :OH => "-", :OD => "--",
                      # "water group" (roughly, I ain't a chemist)
                      :H2O => "-", :HDO => "--",
                      :H2O2 => "-", :HDO2 => "--",
                      :HO2 => "-", :DO2 => "--",
                      # O group
                      :O1D => "-", :O => "-", :O2 => "-", :O3 => "-",
                      # CO group
                      :CO => "-", :CO2 => "-", :CO2pl => "-",
                      :HOCO => "-", :DOCO => "--",
                      # nonreactants
                      :Ar => "-", :N2 => "-",);

################################################################################
##################### FUNDAMENTAL CONSTANTS AND SPECIES LIST ###################
################################################################################

# fundamental constants
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg

# mars parameters
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm

# array of symbols for each species
const fullspecieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
                         :O1D, :H, :N2, :Ar, :CO2pl, :HOCO,
                         # species for deuterium chemistry:
                         :HDO, :OD, :HDO2, :D, :DO2, :HD, :DOCO];
specieslist=fullspecieslist;
const nochemspecies = [:N2, :Ar, :CO2pl, :H2O, :HDO];
const chemspecies = setdiff(specieslist, nochemspecies);
const notransportspecies = [:CO2pl, :H2O, :HDO];
const transportspecies = setdiff(specieslist, notransportspecies);
const speciesmolmasslist = Dict(:CO2=>44, :O2=>32, :O3=>48, :H2=>2, :OH=>17,
                                :HO2=>33, :H2O=>18, :H2O2=>34, :O=>16, :CO=>28,
                                :O1D=>16, :H=>1, :N2=>28, :Ar=>40, :CO2pl=>44,
                                :HOCO=>45,
                                # deuterium chemistry:
                                :HDO=>19, :OD=>18, :HDO2=>35, :D=>2, :DO2=>34,
                                :HD=>3, :DOCO=>46);
                                # we use 40Ar because its the most stable isotope
const fluxlist = map(fluxsymbol, fullspecieslist)

# array of species for which photolysis is important. All rates should
# start with J and contain a species in specieslist above, which is used to 
# compute photolysis. 
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

# Altitude grid discretization
  const zmin = alt[1]
const zmax = alt[end];
const dz = alt[2]-alt[1];

################################################################################
############################### REACTION NETWORK ###############################
################################################################################

# function to replace three body rates with the recommended expression
threebody(k0, kinf) = :($k0*M/(1+$k0*M/$kinf)*0.6^((1+(log10($k0*M/$kinf))^2)^-1))
threebodyca(k0, kinf) = :($k0/(1+$k0/($kinf/M))*0.6^((1+(log10($k0/($kinf*M)))^2)^-1))

# reactions and multipliers on base rates for deuterium reactions from Yung
# 1988; base rates from this work or Chaffin+ 2017. Note: below, H-ana means 
# the same reaction but with only H-bearing species.
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
             [[:O1D, :HDO], [:OD, :OH], :(1.63e-10*exp(60 ./ T))], # Yung88: rate same as H-ana.

             # loss of H2
             [[:H2, :O], [:OH, :H], :(6.34e-12*exp(-4000 ./ T))], # Mike's rate?
             [[:HD, :O], [:OH, :D], :(4.40e-12*exp(-4390 ./ T))], # NIST
             [[:HD, :O], [:OD, :H], :(1.68e-12*exp(-4400 ./ T))], # NIST
             # HD and H2 exchange
             [[:H, :HD], [:H2, :D], :(6.31e-11*exp(-4038 ./ T))], # rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11*exp(-3821 ./ T))], # NIST (1986, 200-300K): 8.19e-13*exp(-2700/T)

             ## OH + H2
             [[:OH, :H2], [:H2O, :H], :(2.8e-12*exp(-1800 ./ T))], # Sander2011. Yung89: 5.5e-12*exp(-2000/T). KIDA: 7.7E-12*exp(-2100/T). old rate from Mike: 9.01e-13*exp(-1526/T)
             [[:OH, :HD], [:HDO, :H], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ T))], # Yung88: rate (3/20)*H-ana. Sander2011: 5e-12*exp(-2130 ./ T)
             [[:OH, :HD], [:H2O, :D], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ T))], # see prev line
             [[:OD, :H2], [:HDO, :H], :(2.8e-12*exp(-1800 ./ T))], # Yung88: rate same as H-ana (assumed)
             [[:OD, :H2], [:H2O, :D], :(0)], # Yung88 (assumed)
             ### [[:OD, :HD], [:HDO, :D], :(???)],  # possibilities for which I 
             ### [[:OD, :HD], [:D2O, :H], :(???)],  # can't find a rate...?

             # recombination of H. Use EITHER the first line OR the 2nd and 3rd.
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32*(298 ./ T)^2.27)],
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
             [[:H, :DO2], [:HO2, :D], :(1e-10/(0.54*exp(890 ./ T)))], # Yung88 (assumed) - turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 - turn off for Case 2

             ## H + H2O2. deuterated analogues added 3/29
             [[:H, :H2O2], [:HO2, :H2],:(2.81e-12*exp(-1890 ./ T))], # verified NIST 4/3/18. Only valid for T>300K. No experiment for lower.
             # [[:H, :HDO2], [:DO2, :H2], :(0)], # Cazaux2010: branching ratio = 0
             # [[:H, :HDO2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:DO2, :H2], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             [[:H, :H2O2], [:H2O, :OH],:(1.7e-11*exp(-1800 ./ T))], # verified NIST 4/3/18
             [[:H, :HDO2], [:HDO, :OH], :(0.5*1.16E-11*exp(-2110 ./ T))], # Cazaux2010: BR = 0.5. Rate for D + H2O2, valid 294<T<464K, NIST, 4/3/18
             [[:H, :HDO2], [:H2O, :OD], :(0.5*1.16E-11*exp(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:HDO, :OH], :(0.5*1.16E-11*exp(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:H2O, :OD], :(0.5*1.16E-11*exp(-2110 ./ T))], # see previous line
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
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5)*8.7e-10)]
             ];

################################################################################
############################ ADD DEUTERATED SPECIES ############################
################################################################################

# D/H ratio ====================================================================
# General D/H ratio for mars, 5.5*SMOW, Atmosphere & Climate of Mars 2017 
if args[1] != "dh"
    DH = 5.5 * 1.6e-4
elseif args[1] == "dh"
    DH = args[2] * 1.6e-4
end

# modify n_current with deuterated species profiles ============================
n_current[:HDO] = n_current[:H2O] * DH
n_current[:OD] = n_current[:OH] * DH
n_current[:HDO2] = n_current[:H2O2] * DH
n_current[:D] = n_current[:H] * DH
n_current[:DO2] = n_current[:HO2] * DH
n_current[:HD] = n_current[:H2] * DH
n_current[:DOCO] = n_current[:HOCO] * DH

# add the new Jrates --the values will get populated automatically =============
n_current[:JHDOtoHpOD] = zeros(length(alt))
n_current[:JHDOtoDpOH] = zeros(length(alt))
n_current[:JHDO2toOHpOD] = zeros(length(alt))
n_current[:JHDOtoHDpO1D] = zeros(length(alt)) 
n_current[:JHDOtoHpDpO] = zeros(length(alt))
n_current[:JODtoOpD] = zeros(length(alt))
n_current[:JHDtoHpD] = zeros(length(alt))
n_current[:JDO2toODpO] = zeros(length(alt))
n_current[:JHDO2toDO2pH] = zeros(length(alt))
n_current[:JHDO2toHO2pD] = zeros(length(alt))
n_current[:JHDO2toHDOpO1D] =  zeros(length(alt))
n_current[:JODtoO1DpD]  =  zeros(length(alt))

################################################################################
####################### TEMPERATURE/PRESSURE PROFILES ##########################
################################################################################

function Tpiecewise(z::Float64, Tsurf, Ttropo, Texo, E="")
    #= DO NOT MODIFY! If you want to change the temperature, define a
    new function or select different arguments and pass to Temp(z)

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 "half-Gaussian" function for temperatures 
    altitudes above the tropopause, with a constant lapse rate (1.4K/km) 
    in the lower atmosphere. The tropopause width is allowed to vary
    in certain cases.

    z: altitude above surface in cm
    Tsurf: Surface temperature in K
    Tropo: tropopause tempearture
    Texo: exobase temperature
    E: type of experiment, used for determining if mesopause width will vary 
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of tropopause. It varies unless we're only varying the 
    # exobase temperature.
    if (E=="tropo") || (E=="surf")
        ztropo_bot = (Ttropo-Tsurf)/(lapserate)
        ztropowidth = ztropo - ztropo_bot
    else
        ztropo_bot = (Ttropo-Tsurf)/(lapserate)
        ztropowidth = ztropo - ztropo_bot
    end

    if z >= ztropo  # upper atmosphere
        return Texo - (Texo - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif ztropo > z >= ztropo - ztropowidth  # tropopause
        return Ttropo
    elseif ztropo-ztropowidth > z  # lower atmosphere
        return Tsurf + lapserate*z
    end
end

# function get_lapserate(Tsurf, Ttropo, lapserate=-1.4e-5, ztropo=90e5)
#     #=
#     Same as Tpiecewise but returns the lapse rate for logging purposes. No,
#     there isn't a better way to do this without it being a PITA, I checked. :(
#     =#
#     # allow varying troposphere width
#     ztropowidth = ztropo - (1/lapserate)*(Ttropo-Tsurf)
#     if ztropowidth < 0   # in case a profile is nonphysical, let lapse rate vary
#         ztropowidth = 30e5
#         m = (ztropo/1e5 - ztropowidth/1e5) / (Ttropo - Tsurf)
#         lapserate = (1/m)*1e-5
#     end
#     return lapserate
# end

function plot_temp_prof(savepath)
    #=
    Input:
        savepath: where to save the temperature profile
    Output: 
        A single panel plot of the temperature profile
    =#

    alt = (0:2e5:zmax)

    fig, ax = subplots(figsize=(4,6))
    ax.set_facecolor("#ededed")
    grid(zorder=0, color="white", which="both")
    for side in ["top", "bottom", "left", "right"]
        ax.spines[side].set_visible(false)
    end

    plot([Temp(a) for a in alt], alt/1e5)

    ax.set_ylabel("Altitude (km)")
    ax.set_yticks([0, 100, zmax/1e5])
    ax.set_yticklabels([0,100,zmax/1e5])
    ax.set_xlabel("Temperature (K)")

    ax.text(Temp(zmax)*0.9, 185, L"T_{exo}")
    ax.text(Temp(100e5)+5, 75, L"T_{tropo}")
    ax.text(Temp(0.0), 10, L"T_{surface}")

    savefig(experimentdir*savepath*"/temp_profile.png", bbox_inches="tight")
end

#If changes to the temperature are needed, they should be made here 
if args[1]=="temp"
    Temp(z::Float64) = Tpiecewise(z, args[2], args[3], args[4])
    Temp_keepSVP(z::Float64) = Tpiecewise(z, 190.0, 110.0, 200.0) # for testing temp without changing SVP.
    # lapserate_logme = get_lapserate(args[2], args[3])
else
    Temp(z::Float64) = Tpiecewise(z, 190.0, 110.0, 200.0)
    # lapserate_logme = get_lapserate(190.0, 110.0)
end

plot_temp_prof(FNext)


################################################################################
############################### WATER PROFILES #################################
################################################################################

# Control functions ============================================================
# 1st term is a conversion factor to convert to (#/cm^3) bceause the 2nd
# term (from Washburn 1924) gives the value in mmHg
Psat(T::Float64) = ((133.3*1e-6)/(boltzmannK * T))*(10^(-2445.5646/T + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 - 6.757169))

function LHDO(T)
    #=
    Latent heat of vaporization of HDO as a function of temperature in K.
    Using Jancso 1974 (https://pubs.acs.org/doi/pdf/10.1021/cr60292a004)
    and Jacobson's Fundamentals of Atmospheric Modeling.

    How I determined this analytical function:
    1. Used Jancso 1974, pg. 734, to get data on the difference of latent heat 
       of evaporation (L_e) for HDO and H2O from 0 to 100°C. 
    2. Used equation 2.54 in Jacobson to get analytical values for L_e of H2O in 
       the same range (0-100°C)
    3. Combined the L_e of H2O with data from Jancso 1974 to get L_HDO.
    4. Noting that L_e expression is roughly linear, fit a line to the L_HDO
       values and extrapolated to Mars temperatures of interest (100 to 350K)

    The data from Jancso 1974 are in cal/mol vs Celsius. We convert to kJ/mol 
    below. 1 kJ = 239 cal => 0.004184 kJ/1 cal. Fit was done in Python in a 
    separate script. parameters below are the output from said fit.
    =#
      
    M = -11.368481383727467  
    a = -2.1942156994200985
    B = 11007.47142671782 # this is cal/mol
    return (M*(T-a) +B) * (0.004184 / 1) # returns in kJ/mol
end

function Psat_HDO(T)
    #=
    Analytical expression for saturation vapor pressure of HDO, using analytical
    latent heat for HDO determined from fitting to data from Jancso 1974.
    The boiling temperature for HDO (100.7°C, 373.85K) and standard pressure are
    used as reference temp/pressure. The gas constant for HDO was calculated as
    R_HDO = kB/m, which probably came from Pierrehumbert's text.

    returns the pressure in #/cm^3. The conversion from N/m^2 to #/cm^3 is 
    divide by N*m (that is, kT) and multiply by the conversion to cm from m 
    (1e-6). 

    Input
        T: a single temperature in K
    Output:
        Pressure in #/cm^3
    =#
    R_HDO = 434.8 # J/kgK
    L_HDO = LHDO(T) * (1000/1) *(1/0.019) # kJ/mol * J/kJ * mol / (kg of HDO)= J / kg
    T0 = 373.85 # boiling temp for liquid HDO
    P0 = 101325 # standard pressure
    Psat_pascals = P0*exp(-(L_HDO/R_HDO)*(1/T - 1/T0))
    # To convert to #/cm^3 from Pa:
    # Pa = Nm^-2 (1/(Nm)) * (1e-6 m^3 / 1 cm^3)
    return Psat_pascals*(1e-6)/(boltzmannK*T) # return in #/cm^3
end

# H2O Water Profile ============================================================
# 
if fix_SVP==true
    H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # for holding SVP fixed
else
    H2Osat = map(x->Psat(x), map(Temp, alt)) # array in #/cm^3 by altitude
end
H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP as fraction of total atmo
# set H2O SVP fraction to minimum for all alts above first time min is reached
H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
               fill(minimum(H2Osatfrac), length(alt)-2-length(H2Oinitfrac))]

# make profile constant in the lower atmosphere (well-mixed).
# when doing water experiments, the temperature profile is such that the minimum
# in the SVP curve occurs at about 55 km alt. This was found by manual tweaking.
# thus we need to only set the lower atmo mixing ratio below that point, or 
# there will be a little spike in the water profile.
if args[1] == "water"
    H2Oinitfrac[findall(x->x<50e5, alt)] .= args[2]
    MR = args[2] # mixing ratio 
else  # i believe 30 km is supposed to be approximately the hygropause.
    H2Oinitfrac[findall(x->x<30e5, alt)] .= 8.1e-4 # 10 pr μm
    MR = 8.1e-4
end

for i in [1:length(H2Oinitfrac);]
    H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
end

# HDO water profile ============================================================
if fix_SVP==true
    HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))  # for holding SVP fixed
else
    HDOsat = map(x->Psat_HDO(x), map(Temp, alt))
end
HDOsatfrac = HDOsat./map(z->n_tot(n_current, z), alt)
# use D/H ratio to set population of HDO
HDOinitfrac = H2Oinitfrac * DH  # initial profile for HDO

# Compute total water column in pr μm starting with mixing ratio array =========
# H2O #/cm^3 = sum(MR * n_tot) for each alt
H2O_per_cc = sum([MR; H2Oinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
HDO_per_cc = sum([MR*DH; HDOinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))

# pr μm = (H2O #/cm^3) * cm * (mol/#) * (H2O g/mol) * (1 cm^3/g) * (10^4 μm/cm)
# where the lone cm is a slice of the atmosphere of thickness dz, #/mol=6.023e23, 
# H2O or HDO g/mol = 18 or 19, cm^3/g = 1 or 19/18 for H2O or HDO.
# written as conversion factors for clarity.
H2Oprum = (H2O_per_cc * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
HDOprum = (HDO_per_cc * dz) * (19/1) * (1/6.02e23) * (19/18) * (1e4/1)

# Detached parcel - water enhancement ==========================================
# add in a detached water vapor parcel, which looks kinda like a gaussian packet
# floating at 60km (Maltagliati 2013). Only for simulations exploring perturbation.
parcel = 1e-6*map(x->80 .* exp(-((x-60)/12.5)^2),alt[2:end-1]/1e5)+H2Oinitfrac
parcel_HDO = parcel * DH

H2O_per_cc_boost = sum([MR; parcel] .* map(z->n_tot(n_current, z), alt[1:end-1]))
HDO_per_cc_boost = sum([MR*DH; parcel_HDO] .* map(z->n_tot(n_current, z), alt[1:end-1]))

H2Oprum_boost = (H2O_per_cc_boost * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
HDOprum_boost = (HDO_per_cc * dz) * (19/1) * (1/6.02e23) * (19/18) * (1e4/1)

# Write out water content to a file ============================================
f = open(experimentdir*FNext*"/water_column_"*FNext*".txt", "w")
write(f, "Total H2O col: $(H2O_per_cc*2e5)\n")
write(f, "Total HDO col: $(HDO_per_cc*2e5)\n")
write(f, "Total water col: $((H2O_per_cc + HDO_per_cc)*2e5)\n")
write(f, "H2O+HDO at surface: $((H2O_per_cc[1] + HDO_per_cc[1])*2e5)\n")
write(f, "Total H2O (pr μm): $(H2Oprum)\n")
write(f, "Total H2O + parcel (pr μm): $(H2Oprum_boost) \n")
write(f, "Total HDO (pr μm): $(HDOprum)\n")
write(f, "Total HDO + parcel (pr μm): $(HDOprum_boost) \n")
write(f, "Total H2O+HDO, no enhancement: $(H2Oprum + HDOprum)")
close(f)

# Plot the water profiles (as raw mixing ratio) ================================
fig, ax = subplots(figsize=(6,9))
ax.set_facecolor("#ededed")
grid(zorder=0, color="white", which="both")
for side in ["top", "bottom", "left", "right"]
    ax.spines[side].set_visible(false)
end
semilogx(H2Oinitfrac, alt[2:end-1]/1e5, color="cornflowerblue", linewidth=3, 
         label=L"H$_2$O")
semilogx(HDOinitfrac, alt[2:end-1]/1e5, color="cornflowerblue", linestyle="--", 
         linewidth=3, label="HDO")
semilogx(H2Osatfrac, alt[1:end]/1e5, color="black", alpha=0.5, linewidth=3, 
         label=L"H$_2$O saturation")
semilogx(HDOsatfrac, alt[1:end]/1e5, color="black", alpha=0.5, linestyle="--", 
         linewidth=3, label="HDO saturation")
xlabel("Mixing ratio", fontsize=18)
ylabel("Altitude [km]", fontsize=18)
title(L"H$_2$O and HDO model profiles", fontsize=20)
ax.tick_params("both",labelsize=16)
legend()
savefig(experimentdir*FNext*"/water_profile_rawMR.png")

################################################################################
############################# BOUNDARY CONDITIONS ##############################
################################################################################

H_effusion_velocity = effusion_velocity(Temp(zmax), 1.0)
H2_effusion_velocity = effusion_velocity(Temp(zmax), 2.0)
D_effusion_velocity = effusion_velocity(Temp(zmax), 2.0)
HD_effusion_velocity = effusion_velocity(Temp(zmax), 3.0)

#=
    boundary conditions for each species (mostly from Nair 1994, Yung 1988). For 
    most species, default boundary condition is zero (f)lux at top and bottom. 
    Atomic/molecular hydrogen and deuterated analogues have a nonzero effusion 
    (v)elocity at the upper layer of the atmosphere. Some species have a (n)umber 
    density boundary condition.
=#

if args[1] == "Oflux"
    Of = args[2]
else
    Of = 1.2e8
end

const speciesbclist=Dict(
                :CO2=>["n" 2.1e17; "f" 0.],
                :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                :HDO=>["n" HDOsat[1]; "f" 0.],
                :O=>["f" 0.; "f" Of],
                :H2=>["f" 0.; "v" H2_effusion_velocity],
                :HD=>["f" 0.; "v" HD_effusion_velocity],
                :H=>["f" 0.; "v" H_effusion_velocity],
                :D=>["f" 0.; "v" D_effusion_velocity],
               );

################################################################################
############################ DIFFUSION COEFFICIENTS ############################
################################################################################

#=
    molecular diffusion parameters. value[1] = A, value[2] = s in the equation
    D = AT^s / n given by Banks & Kockarts Aeronomy part B eqn. 15.30 and Hunten
    1973, Table 1.

    molecular diffusion is different only for small molecules and atoms
    (H, D, HD, and H2), otherwise all species share the same values (Krasnopolsky
    1993 <- Hunten 1973; Kras cites Banks & Kockarts, but this is actually an
    incorrect citation.)

    D and HD params are estimated by using Banks & Kockarts eqn 15.29 (the
    coefficient on T^0.5/n) to calculate A for H and H2 and D and HD, then using
    (A_D/A_H)_{b+k} = (A_D/A_H)_{hunten} to determine (A_D)_{hunten} since Hunten
    1973 values are HALF what the calculation provides.

    s, the power of T, is not calculated because there is no instruction on how to
    do so and is assumed the same as for the hydrogenated species.

    THESE ARE IN cm^-2 s^-2!!!
=#

diffparams(species) = get(Dict(:H=>[8.4, 0.597], :H2=>[2.23, 0.75],
                               :D=>[5.98, 0.597], :HD=>[1.84, 0.75]),
                               species,[1.0, 0.75])

# override to use altitude instead of temperature
Dcoef(z, species::Symbol, n_current) = Dcoef(Temp(z),n_tot(n_current, z),species)

# thermal diffusion factors (from Krasnopolsky 2002)
# TODO: these probably need to be verified for deuterated species.
thermaldiff(species) = get(Dict(:H=>-0.25, :H2=>-0.25, :D=>-0.25, :HD=>-0.25,
                                :He=>-0.25), species, 0)


################################################################################
####################### COMBINED CHEMISTRY AND TRANSPORT #######################
################################################################################

#=
    We now have objects that return the list of indices and coefficients
    for transport, assuming no other species in the atmosphere
    (transportmat), and for chemistry, assuming no other altitudes
    (chemical_jacobian). We need to perform a kind of outer product on
    these operators, to determine a fully coupled set of equations for
    all species at all altitudes.
=#

# need to get a list of all species at all altitudes to iterate over
const intaltgrid = round.(Int64, alt/1e5)[2:end-1]
const replacespecies = [fullspecieslist, Jratelist, [:T,:M]]

#=
    the rates at each altitude can be computed using the reaction network
    already in place, plus additional equations describing the transport
    to and from the cells above and below:
=#
upeqns = [Any[Any[[s], [Symbol(string(s)*"_above")],Symbol("t"*string(s)*"_up")],
          Any[[Symbol(string(s)*"_above")],[s],Symbol("t"*string(s)*"_above_down")]]
          for s in specieslist]

downeqns = [Any[Any[[s], [Symbol(string(s)*"_below")],Symbol("t"*string(s)*"_down")],
            Any[[Symbol(string(s)*"_below")],[s],Symbol("t"*string(s)*"_below_up")]]
            for s in specieslist]

local_transport_rates = [[[Symbol("t"*string(s)*"_up") for s in specieslist]
                          [Symbol("t"*string(s)*"_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_above_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_below_up") for s in specieslist]]...;]

transportnet = [[upeqns...;]; [downeqns...;]]

# define names for all the species active in the coupled rates:
activespecies = union(chemspecies, transportspecies)
active_above = [Symbol(string(s)*"_above") for s in activespecies]
active_below = [Symbol(string(s)*"_below") for s in activespecies]
inactivespecies = intersect(nochemspecies, notransportspecies)


# obtain the rates and jacobian for each altitude
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x),activespecies)...);
const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

arglist_local = [activespecies; active_above; active_below; inactivespecies;
                 Jratelist; :T; :M; local_transport_rates; :dt]

arglist_local_typed=[:($s::Float64) for s in arglist_local]

@eval begin
    function ratefn_local($(arglist_local_typed[1:end-1]...))
        $rates_local # evaluates the rates_local expression
    end
end

@eval begin
    function chemJmat_local($(arglist_local_typed...))
        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])
        localchemJval = -dt*$(chemJ_local[3])

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = -dt*$(chemJ_above[3])

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = -dt*$(chemJ_below[3])

        ((localchemJi, localchemJj, localchemJval),
         (abovechemJi, abovechemJj, abovechemJval),
         (belowchemJi, belowchemJj, belowchemJval))
    end
end

@eval begin
    function reactionrates_local($(specieslist...), $(Jratelist...), T, M)
        #= a function to return chemical reaction rates for specified species
           concentrations =#
        $(Expr(:vcat, map(x->Expr(:call,:*,x[1]..., x[3]), reactionnet)...))
    end
end

################################################################################
######################### PHOTOCHEMICAL CROSS SECTIONS #########################
################################################################################

# Change following line as needed depending on local machine
xsecfolder = scriptdir * "uvxsect/";

# Crosssection Files ===========================================================
co2file = "CO2.dat"
co2exfile = "binnedCO2e.csv"
h2ofile = "h2oavgtbl.dat"
hdofile = "HDO.dat"
h2o2file = "H2O2.dat"
hdo2file = "H2O2.dat" #TODO: do HDO2 xsects exist?
o3file = "O3.dat"
o3chapfile = "O3Chap.dat"
o2file = "O2.dat"
o2_130_190 = "130-190.cf4"
o2_190_280 = "190-280.cf4"
o2_280_500 = "280-500.cf4"
h2file = "binnedH2.csv"
hdfile = "binnedH2.csv" # TODO: change this to HD file if xsects ever exist
ohfile = "binnedOH.csv"
oho1dfile = "binnedOHo1D.csv"
odfile = "OD.csv"

# Loading Data =================================================================
# CO2 --------------------------------------------------------------------------
# temperature-dependent between 195-295K
co2xdata = readandskip(xsecfolder*co2file,'\t',Float64, skipstart = 4)
function co2xsect(T::Float64)
    clamp(T, 195, 295)
    Tfrac = (T-195)/(295-195)

    arr = [co2xdata[:,1]; (1-Tfrac)*co2xdata[:,2]+Tfrac*co2xdata[:,3]]
    reshape(arr, length(co2xdata[:,1]),2)
end

# CO2 photoionization (used to screen high energy sunlight)
co2exdata = readandskip(xsecfolder*co2exfile,',',Float64, skipstart = 4)

# H2O & HDO --------------------------------------------------------------------
h2oxdata = readandskip(xsecfolder*h2ofile,'\t',Float64, skipstart = 4)

# These crosssections for HDO are for 298K.
hdoxdata = readandskip(xsecfolder*hdofile,'\t', Float64, skipstart=12)

# H2O2 + HDO2 ------------------------------------------------------------------
# from 260-350 the following analytic calculation fitting the
# temperature dependence is recommended by Sander 2011:
function h2o2xsect_l(l::Float64, T::Float64)
    #=
    Analytic calculation of H2O2 cross section using temperature dependencies
    l: wavelength in nm
    T: temperature in K
    =#
    l = clamp(l, 260, 350)
    T = clamp(T, 200, 400)

    A = [64761., -921.70972, 4.535649,
       -0.0044589016, -0.00004035101,
       1.6878206e-7, -2.652014e-10, 1.5534675e-13]
    B = [6812.3, -51.351, 0.11522, -0.000030493, -1.0924e-7]

    lpowA = map(n->l^n,[0:7;])
    lpowB = map(n->l^n,[0:4;])

    expfac = 1.0/(1+exp(-1265/T))

    return 1e-21*(expfac*reduce(+, map(*, A, lpowA))+(1-expfac)*reduce(+, map(*, B, lpowB)))
end

function h2o2xsect(T::Float64)
    #=
    stitches together H2O2 cross sections, some from Sander 2011 table and some
    from the analytical calculation recommended for 260-350nm recommended by the
    same.
    T: temperature in K
    =#
    retl = h2o2xdata[:,1]
    retx = 1e4*h2o2xdata[:,2]#factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    return reshape([retl; retx], length(retl), 2)
end

function hdo2xsect(T::Float64)
    retl = hdo2xdata[:,1]
    retx = 1e4*hdo2xdata[:,2] #factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    reshape([retl; retx], length(retl), 2)
end

# the data in the following table cover the range 190-260nm
h2o2xdata = readandskip(xsecfolder*h2o2file,'\t',Float64, skipstart=3)
hdo2xdata = readandskip(xsecfolder*hdo2file,'\t',Float64, skipstart=3)

# O3 ---------------------------------------------------------------------------
# including IR bands which must be resampled from wavenumber
o3xdata = readandskip(xsecfolder*o3file,'\t',Float64, skipstart=3)
o3ls = o3xdata[:,1]
o3xs = o3xdata[:,2]
o3chapxdata = readandskip(xsecfolder*o3chapfile,'\t',Float64, skipstart=3)
o3chapxdata[:,1] = map(p->1e7/p, o3chapxdata[:,1])
for i in [round(Int, floor(minimum(o3chapxdata[:,1]))):round(Int, ceil(maximum(o3chapxdata))-1);]
    posss = getpos(o3chapxdata, x->i<x<i+1)
    dl = diff([map(x->o3chapxdata[x[1],1], posss); i])
    x = map(x->o3chapxdata[x[1],2],posss)
    ax = reduce(+,map(*,x, dl))/reduce(+,dl)
    global o3ls = [o3ls; i+0.5]
    global o3xs = [o3xs; ax]
end
o3xdata = reshape([o3ls; o3xs],length(o3ls),2)

# O2 ---------------------------------------------------------------------------
# including temperature-dependent Schumann-Runge bands.
function binupO2(list)
    ret = Float64[];
    for i in [176:203;]
        posss = getpos(list[:,1],x->i<x<i+1)
        dl = diff([map(x->list[x[1],1],posss); i])
        x0 = map(x->list[x[1],2],posss)
        x1 = map(x->list[x[1],3],posss)
        x2 = map(x->list[x[1],4],posss)
        ax0 = reduce(+,map(*,x0, dl))/reduce(+,dl)
        ax1 = reduce(+,map(*,x1, dl))/reduce(+,dl)
        ax2 = reduce(+,map(*,x2, dl))/reduce(+,dl)
        append!(ret,[i+0.5, ax0, ax1, ax2])
    end
    return transpose(reshape(ret, 4, 203-176+1))
end

function o2xsect(T::Float64)
    o2x = deepcopy(o2xdata);
    #fill in the schumann-runge bands according to Minschwaner 1992
    T = clamp(T, 130, 500)
    if 130<=T<190
        o2schr = o2schr130K
    elseif 190<=T<280
        o2schr = o2schr190K
    else
        o2schr = o2schr280K
    end

    del = ((T-100)/10)^2

    for i in [176.5:203.5;]
        posO2 = something(findfirst(isequal(i), o2x[:, 1]), 0)
        posschr = something(findfirst(isequal(i), o2schr[:, 1]), 0)
        o2x[posO2, 2] += 1e-20*(o2schr[posschr, 2]*del^2
                                + o2schr[posschr, 3]*del
                                + o2schr[posschr, 4])
    end

    # add in the herzberg continuum (though tiny)
    # measured by yoshino 1992
    for l in [192.5:244.5;]
        posO2 = something(findfirst(isequal(l), o2x[:, 1]), 0)
        o2x[posO2, 2] += 1e-24*(-2.3837947e4
                            +4.1973085e2*l
                            -2.7640139e0*l^2
                            +8.0723193e-3*l^3
                            -8.8255447e-6*l^4)
    end

    return o2x
end

o2xdata = readandskip(xsecfolder*o2file,'\t',Float64, skipstart = 3)
o2schr130K = readandskip(xsecfolder*o2_130_190,'\t',Float64, skipstart = 3)
o2schr130K[:,1] = map(p->1e7/p, o2schr130K[:,1])
o2schr130K = binupO2(o2schr130K)
o2schr190K = readandskip(xsecfolder*o2_190_280,'\t',Float64, skipstart = 3)
o2schr190K[:,1] = map(p->1e7/p, o2schr190K[:,1])
o2schr190K = binupO2(o2schr190K)
o2schr280K = readandskip(xsecfolder*o2_280_500,'\t',Float64, skipstart = 3)
o2schr280K[:,1] = map(p->1e7/p, o2schr280K[:,1])
o2schr280K = binupO2(o2schr280K)

# HO2 & DO2 --------------------------------------------------------------------
function ho2xsect_l(l::Float64)
    #= function to compute HO2 cross-section as a function of wavelength l
    in nm, as given by Sander 2011 JPL Compilation =#
    a = 4.91
    b = 30612.0
    sigmamed = 1.64e-18
    vmed = 50260.0
    v = 1e7/l;
    if 190<=l<=250
        return HO2absx = sigmamed / ( 1 - b/v ) * exp( -a * log( (v-b)/(vmed-b) )^2 )
    else
        return 0.0
    end
end

ho2xsect = [190.5:249.5;]
ho2xsect = reshape([ho2xsect; map(ho2xsect_l, ho2xsect)],length(ho2xsect),2)
do2xsect = deepcopy(ho2xsect) # TODO: find crosssection for DO2

# H2 & HD ----------------------------------------------------------------------
h2xdata = readandskip(xsecfolder*h2file,',',Float64, skipstart=4)
hdxdata = readandskip(xsecfolder*hdfile,',',Float64, skipstart=4)

# OH & OD ----------------------------------------------------------------------
ohxdata = readandskip(xsecfolder*ohfile,',',Float64, skipstart=4)
ohO1Dxdata = readandskip(xsecfolder*oho1dfile,',',Float64, skipstart=4)
odxdata = readandskip(xsecfolder*odfile,',',Float64, skipstart=3)

# PHOTODISSOCIATION ============================================================
function padtosolar(crosssection::Array{Float64, 2})
    # a function to take an Nx2 array and pad it with zeroes until it's the
    # same length as the solar flux. Returns the cross sections only, as
    # the wavelengths are shared by solarflux
    positions = map(x->something(findfirst(isequal(x), solarflux[:,1]), 0), crosssection[:,1])
    retxsec = fill(0.,length(solarflux[:,1]))
    retxsec[positions] = crosssection[:,2]
    return retxsec
end

function quantumyield(xsect::Array, arr)
    #= 
        function to assemble cross-sections for a given pathway. Inputs are
        an Nx2 array xsect with wavelengths and photoabsorption cross
        sections, and arr, a tuple of tuples with a condition and a quantum
        yield multiplicative factor, either constant or a function of
        wavelength in the given regime. Return is an array with all of the
        matching wavelengths and the scaled cross-sections.
    =#
    lambdas = Float64[];
    rxs = Float64[];
    for (cond, qeff) in arr
        places = findall(cond, xsect[:,1])
        append!(lambdas, xsect[places, 1])
        #if we have a number then map to a function
        isa(qeff, Function) ? (qefffn = qeff) : (qefffn = x->qeff)
        append!(rxs, map(*,map(qefffn, xsect[places, 1]),xsect[places, 2]))
    end

    return reshape([lambdas; rxs],length(lambdas),2)
end

# Change following line as needed depending on local machine
const solarflux=readandskip(scriptdir*"marssolarphotonflux.dat",'\t', 
                            Float64,skipstart=4)[1:2000,:]
solarflux[:,2] = solarflux[:,2]/2

absorber = Dict(:JCO2ion =>:CO2,
                :JCO2toCOpO =>:CO2,
                :JCO2toCOpO1D =>:CO2,
                :JO2toOpO =>:O2,
                :JO2toOpO1D =>:O2,
                :JO3toO2pO =>:O3,
                :JO3toO2pO1D =>:O3,
                :JO3toOpOpO =>:O3,
                :JH2toHpH =>:H2,
                :JHDtoHpD => :HD,
                :JOHtoOpH =>:OH,
                :JOHtoO1DpH =>:OH,
                :JODtoOpD =>:OD,
                :JODtoO1DpD => :OD,
                :JHO2toOHpO =>:HO2,
                :JDO2toODpO => :DO2,
                :JH2OtoHpOH =>:H2O,
                :JH2OtoH2pO1D =>:H2O,
                :JH2OtoHpHpO =>:H2O,
                :JH2O2to2OH =>:H2O2,
                :JH2O2toHO2pH =>:H2O2,
                :JH2O2toH2OpO1D =>:H2O2,
                :JHDO2toHDOpO1D => :HDO2,
                :JHDOtoHpOD=>:HDO,
                :JHDO2toOHpOD=>:HDO2,
                :JHDO2toDO2pH => :HDO2,
                :JHDO2toHO2pD => :HDO2,
                :JHDOtoDpOH=>:HDO,
                :JHDOtoHpDpO=>:HDO,
                :JHDOtoHDpO1D=>:HDO
                );

#=
    this is a dictionary of the 1-nm photodissociation or photoionization
    cross-sections important in the atmosphere. keys are symbols found in
    jratelist. each entry is an array of arrays, yielding the wavelengths
    and cross-sections for each altitude in the atmosphere.

    NOTE: jspecies refers to the photodissociation or photoionization
    cross section for a particular species which produces a UNIQUE SET OF
    PRODUCTS. In this sense, crosssection has already folded in quantum
    efficiency considerations.
=#
crosssection = Dict{Symbol, Array{Array{Float64}}}()
# now add the cross-sections

# CO2 photodissociation --------------------------------------------------------
setindex!(crosssection, fill(co2exdata, length(alt)), :JCO2ion)
#CO2+hv->CO+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->l>167, 1), (l->95>l, 0.5))),
          map(t->co2xsect(t),map(Temp, alt))), :JCO2toCOpO)
#CO2+hv->CO+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->95<l<167, 1), (l->l<95, 0.5))),
          map(t->co2xsect(t),map(Temp, alt))), :JCO2toCOpO1D)

# O2 photodissociation ---------------------------------------------------------
#O2+hv->O+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x>175, 1),)), map(t->o2xsect(t),map(Temp, alt))),
          :JO2toOpO)
#O2+hv->O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<175, 1),)), map(t->o2xsect(t),map(Temp, alt))),
          :JO2toOpO1D)

# O3 photodissociation ---------------------------------------------------------
# The quantum yield of O1D from ozone photolysis is actually
# well-studied! This adds some complications for processing.
function O3O1Dquantumyield(lambda, temp)
    if lambda < 306. || lambda > 328.
        return 0.
    end
    temp=clamp(temp, 200, 320)#expression is only valid in this T range

    X = [304.225, 314.957, 310.737];
    w = [5.576, 6.601, 2.187];
    A = [0.8036, 8.9061, 0.1192];
    v = [0.,825.518];
    c = 0.0765;
    R = 0.695;
    q = exp.(-v/(R*temp))
    qrat = q[1]/(q[1]+q[2])

    (q[1]/sum(q)*A[1]*exp.(-((X[1]-lambda)/w[1])^4.)
     +q[2]/sum(q)*A[2]*(temp/300.)^2 .* exp.(-((X[2]-lambda)/w[2])^2.)
     +A[3]*(temp/300.)^1.5*exp.(-((X[3]-lambda)/w[3])^2.)
     +c)
end

# O3+hv->O2+O
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1-(1.37e-2*193-2.16)),
                               (l->193<=l<225, l->(1 .- (1.37e-2*l-2.16))),
                               (l->225<=l<306, 0.1),
                               (l->306<=l<328, l->(1 .- O3O1Dquantumyield(l, t))),
                               (l->328<=l<340, 0.92),
                               (l->340<=l, 1.0)
                               ))
              ,map(Temp, alt)
              ),
          :JO3toO2pO)
# O3+hv->O2+O1D
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1.37e-2*193-2.16),
                               (l->193<=l<225, l->(1.37e-2*l-2.16)),
                               (l->225<=l<306, 0.9),
                               (l->306<=l<328, l->O3O1Dquantumyield(l, t)),
                               (l->328<=l<340, 0.08),
                               (l->340<=l, 0.0)
                               ))
              ,map(Temp, alt)
              ),
          :JO3toO2pO1D)
# O3+hv->O+O+O
setindex!(crosssection,
          fill(quantumyield(o3xdata,((x->true, 0.),)),length(alt)),
          :JO3toOpOpO)

# H2 and HD photodissociation --------------------------------------------------
# H2+hv->H+H
setindex!(crosssection, fill(h2xdata, length(alt)), :JH2toHpH)
# HD+hν -> H+D  # TODO: Do we need to x JHDtoHpD by some factor?
setindex!(crosssection, fill(hdxdata, length(alt)), :JHDtoHpD)

# OH and OD photodissociation --------------------------------------------------
# OH+hv->O+H
setindex!(crosssection, fill(ohxdata, length(alt)), :JOHtoOpH)
# OH+hv->O1D+H
setindex!(crosssection, fill(ohO1Dxdata, length(alt)), :JOHtoO1DpH)
# OD + hv -> O+D  # TODO: do we need to x JODtoOpD by some factor?
setindex!(crosssection, fill(odxdata, length(alt)), :JODtoOpD)
# OD + hν -> O(¹D) + D  # TODO: do we need to x JODtoO1DpD by some factor?
setindex!(crosssection, fill(ohO1Dxdata, length(alt)), :JODtoO1DpD)

# HO2 and DO2 photodissociation ------------------------------------------------
# HO2 + hν -> OH + O
setindex!(crosssection, fill(ho2xsect, length(alt)), :JHO2toOHpO)
# DO2 + hν -> OD + O   # TODO: do we need to x JDO2toODpO by some factor?
setindex!(crosssection, fill(do2xsect, length(alt)), :JDO2toODpO)

# H2O and HDO photodissociation ------------------------------------------------
# H2O+hv->H+OH
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),length(alt)),
          :JH2OtoHpOH)

# H2O+hv->H2+O1D
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JH2OtoH2pO1D)

# H2O+hv->H+H+O
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->true, 0),)),length(alt)),
          :JH2OtoHpHpO)

# HDO + hν -> H + OD
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alt)),
          :JHDOtoHpOD)

# HDO + hν -> D + OH
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alt)),
          :JHDOtoDpOH)

# HDO + hν -> HD + O1D
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JHDOtoHDpO1D)

# HDO + hν -> H + D + O
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->true, 0),)),length(alt)),
          :JHDOtoHpDpO)


# H2O2 and HDO2 photodissociation ----------------------------------------------
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
          map(t->h2o2xsect(t), map(Temp, alt))), :JH2O2to2OH)

# H2O2+hv->HO2+H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.15),(x->x>230, 0))),
          map(t->h2o2xsect(t), map(Temp, alt))), :JH2O2toHO2pH)

# H2O2+hv->H2O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->h2o2xsect(t),
          map(Temp, alt))), :JH2O2toH2OpO1D)

# HDO2 + hν -> OH + OD
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
          map(t->hdo2xsect(t), map(Temp, alt))), :JHDO2toOHpOD)

# HDO2 + hν-> DO2 + H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
          map(t->hdo2xsect(t), map(Temp, alt))), :JHDO2toDO2pH)

# HDO2 + hν-> HO2 + D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
          map(t->hdo2xsect(t), map(Temp, alt))), :JHDO2toHO2pD)

# HDO2 + hν -> HDO + O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->hdo2xsect(t),
          map(Temp, alt))), :JHDO2toHDOpO1D)

# Solar Input ------------------------------------------------------------------
lambdas = Float64[]
for j in Jratelist, ialt in 1:length(alt)
    global lambdas = union(lambdas, crosssection[j][ialt][:,1])
end

if !(setdiff(solarflux[:,1],lambdas)==[])
    throw("Need a broader range of solar flux values!")
end

# pad all cross-sections to solar
for j in Jratelist, ialt in 1:length(alt)
    crosssection[j][ialt] = padtosolar(crosssection[j][ialt])
end

# we need some global objects for the Jrates calculation:
# intensity as a function of wavelength at each altitude

# this is the unitialized array for storing values
solarabs = fill(fill(0.,size(solarflux, 1)),length(alt)-2);

################################################################################
############################# CONVERGENCE CODE #################################
################################################################################

# set the water profiles =======================================================
n_current[:H2O] = H2Oinitfrac.*map(z->n_tot(n_current, z), alt[2:end-1])
n_current[:HDO] = HDOinitfrac.*map(z->n_tot(n_current, z), alt[2:end-1])

# Plot initial water profile ===================================================
# ...in ppm --------------------------------------------------------------------
fig, ax = subplots(figsize=(6,8))
ax.set_facecolor("#ededed")
grid(zorder=0, color="white", which="both")
for side in ["top", "bottom", "left", "right"]
    ax.spines[side].set_visible(false)
end
ax.tick_params(axis="x", which="minor", bottom=true, top=true)
semilogx(H2Oinitfrac/1e-6, alt[2:end-1]/1e5, color="cornflowerblue", linewidth=2,
         label=L"H$_2$O")
semilogx(HDOinitfrac/1e-6, alt[2:end-1]/1e5, color="cornflowerblue", 
         linestyle="--", linewidth=2, label="HDO")
xlabel("Volume Mixing Ratio [ppm]")
ylabel("Altitude [km]")
title(L"H$_2$O and HDO model profiles")
legend()
savefig(experimentdir*FNext*"/initial_water_MR.png")

# ...and raw abundance ---------------------------------------------------------
fig, ax = subplots(figsize=(6,8))
ax.set_facecolor("#ededed")
grid(zorder=0, color="white", which="both")
for side in ["top", "bottom", "left", "right"]
    ax.spines[side].set_visible(false)
end
ax.tick_params(axis="x", which="minor", bottom=true, top=true)
semilogx(n_current[:H2O], alt[2:end-1]/1e5, color="cornflowerblue", linewidth=2,
         label=L"H$_2$O")
semilogx(n_current[:HDO], alt[2:end-1]/1e5, color="cornflowerblue", 
         linestyle="--", linewidth=2, label="HDO")
xlabel(L"Species density [cm$^{-2}$")
ylabel("Altitude [km]")
title(L"H$_2$O and HDO model profiles")
legend()
savefig(experimentdir*FNext*"/initial_water_number.png")
#show()

# do the convergence ===========================================================
# initialize whole atmosphere figure 
fig, ax = subplots(figsize=(10,6))
ax.set_facecolor("#ededed")
ax.grid(color="white")
for side in ["top", "bottom", "left", "right"]
    ax.spines[side].set_visible(false)
end

[timeupdate(t) for t in [10.0^(1.0*i) for i in -3:14]]
for i in 1:100
    plotatm()
    # println("dt: 1e14 iter $(i)")
    update!(n_current, 1e14)
end

# write out the new converged file to matching folder. 
towrite = experimentdir*FNext*"/converged_standardwater_D_"*FNext*".h5"
write_ncurrent(n_current, towrite)
println("Wrote $(towrite)")

# save the figure
savefig(experimentdir*FNext*"/converged_"*FNext*".png", bbox_inches="tight")
println("Saved figure to same folder")

################################################################################
################################# LOGGING ######################################
################################################################################

# crosssection dict for logging purposes =======================================
xsect_dict = Dict("CO2"=>[co2file, co2exfile], 
              "H2O, HDO"=>[h2ofile, hdofile],
              "H2O2, HDO2"=>[h2o2file, hdo2file],
              "O3"=>[o3file, o3chapfile],
              "O2"=>[o2file, o2_130_190, o2_190_280, o2_280_500],
              "H2, HD"=>[h2file, hdfile],
              "OH, OD"=>[ohfile, oho1dfile, odfile])

# Log temperature and water parameters =========================================
if args[1]=="temp"
    input_string = "T_0=$(args[2]), T_tropo=$(args[3]), T_exo=$(args[4])" * 
                   "\nwater init=8e-4\nDH=5.5 \nOflux=1.2e8" #*
                   #"\nlapse rate=$(lapserate_logme)\n"
elseif args[1]=="water"
    input_string = "T_0=190.0, T_tropo=110.0, T_exo=200.0\n" *
                   "water init=$(args[2])\nDH=5.5\nOflux=1.2e8\n" #*
                   #"lapse rate=$(lapserate_logme)\n"
elseif args[1]=="dh"
    input_string = "T_0=190.0, T_tropo=110.0, T_exo=200.0\nwater=8e-4\n" *
                   "DH=$(args[2]) \nOflux=1.2e8\n"#lapse rate=$(lapserate_logme)\n"
elseif args[1]=="Oflux"
    input_string = "T_0=190.0, T_tropo=110.0, T_exo=200.0 \nwater=8e-4" *
                   "\nDH=5.5\nOflux=$(Of)\n"#lapse rate=$(lapserate_logme)\n"
end

# Write the log ================================================================
f = open(experimentdir*FNext*"/convergence_"*FNext*".txt", "w")
write(f, "Finished convergence for $(args[1]) experiment with control parameters: \n")
write(f, input_string)
write(f, "\nSVP fixed: $(fix_SVP)\n")
write(f, "\nCROSS SECTIONS: \n")
for k in keys(xsect_dict)  # cross sections
    write(f, k*": "*join(xsect_dict[k], ", ")*"\n")
end
# boundary conditions
write(f, "\nBOUNDARY CONDITIONS: \n")
write(f, "n: number density at surface, f: flux at top, v: velocity at top\n")
for k2 in keys(speciesbclist)
    bcstring = join([join(speciesbclist[k2][1, :], "="), 
                     join(speciesbclist[k2][2, :], "=")], ", ")
    write(f, string(k2)*": "*bcstring*"\n")
end

#write(f, "note: this simulation was run with SOLAR MAX values")
# write(f, "note: This simulation converged with SVP for mean temperatures, but with 
#           varying temperature profile")
close(f)

println("ALERT: Finished convergence")
println()

