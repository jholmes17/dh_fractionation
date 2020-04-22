module Photochemistry

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using PlotUtils

include("../../.././PARAMETERS.jl")

export create_folder, getpos, deletefirst, fluxsymbol, readandskip,
       get_ncurrent, write_ncurrent, n_tot, 
       #scaleH, 
       #fluxcoefs, lower_up, upper_down, 
       meanmass, loss_equations,loss_rate, production_equations, production_rate,
           chemical_jacobian, getrate, getflux, fluxes, ratefn, chemJmat, #reactionrates, 
       plotatm, 
       Keddy, Dcoef, 
       effusion_velocity, speciesbcs, boundaryconditions,
       Tpiecewise, Psat, Psat_HDO

# auxiliary functions ==========================================================

function create_folder(foldername, parentdir)
    #=
    Creates a folder, or if it already exists, notifies the user.
    =#
    println("Checking for existence of $(foldername) folder in $(parentdir)")
    dircontents = readdir(parentdir)
    if foldername in dircontents
        println("Folder $(foldername) exists")
    else
        mkdir(parentdir*foldername)
        println("Created folder ", foldername)
    end
end

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

function fluxsymbol(x)
    Symbol(string("f",string(x)))
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


# Main array functions =========================================================

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

function write_ncurrent(n_current, filename, alt)
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

function n_tot(n_current, z)
    #= get the total number density at a given altitude =#
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

function n_tot(n_current, z, n_alt_index)
    #= get the total number , density at a given altitude =#
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end


# transport/scale height =======================================================

# function scaleH(z, T::Float64, mm::Real)
#     #= Computes the scale height of the atmosphere for the mean molar mass =#
#     return boltzmannK*T/(mm*mH*marsM*bigG)*(((z+radiusM)*1e-2)^2)*1e2
#     # constants are in MKS. Convert to m and back to cm.
# end

# function scaleH(z, species::Symbol)
#     #=
#     computes scale height of the atmosphere for a specific species
#     =#
#     T=Temp(z)
#     mm = speciesmolmasslist[species]
#     scaleH(z, T, mm)
# end

# function scaleH(z, T::Float64, species::Symbol)
#     #=
#     computes scale height of the atmosphere for a specific species
#     =#
#     mm = speciesmolmasslist[species]
#     scaleH(z, T, mm)
# end

# function scaleH(z, T::Float64, n_current)
#     #= Computes the scale height of the atmosphere for the mean molar mass =#
#     mm = meanmass(n_current, z)
#     scaleH(z, T, mm)
# end

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

# function fluxcoefs(z, dz, Kv, Dv, Tv, Hsv, H0v, species) #ntv,
#     #= 
#     function to generate coefficients of the transport network.
#     z: altitude
#     dz: altitude layer thickness ("resolution")
#     Kv: eddy diffusion coefficient
#     Dv: molecular diffusion coefficient
#     Tv: temperature
#     Hsv: scale height by species
#     H0v: mean atmospheric scale height
#     species: species symbol 

#     v refers to "vector"
#     u refers to "upper" (the upper parcel)
#     l refers to "lower" (the lower parcel)
#     =#
#     Dl = (Dv[1] + Dv[2])/2.0
#     Kl = (Kv[1] + Kv[2])/2.0
#     Tl = (Tv[1] + Tv[2])/2.0
#     dTdzl = (Tv[2] - Tv[1])/dz
#     Hsl = (Hsv[1] + Hsv[2])/2.0
#     H0l = (H0v[1] + H0v[2])/2.0

#     # we have two flux terms to combine:
#     sumeddyl = (Dl+Kl)/dz/dz
#     gravthermall = (Dl*(1/Hsl + (1+thermaldiff(species))/Tl*dTdzl) +
#                     Kl*(1/H0l + 1/Tl*dTdzl))/(2*dz)

#     Du = (Dv[2] + Dv[3])/2.0
#     Ku = (Kv[2] + Kv[3])/2.0
#     Tu = (Tv[2] + Tv[3])/2.0
#     dTdzu = (Tv[3] - Tv[2])/dz
#     Hsu = (Hsv[2] + Hsv[3])/2.0
#     H0u = (H0v[2] + H0v[3])/2.0

#     # we have two flux terms to combine:
#     sumeddyu = (Du+Ku)/dz/dz
#     gravthermalu = (Du*(1/Hsu + (1 + thermaldiff(species))/Tu*dTdzu) +
#                   Ku*(1/H0u + 1/Tu*dTdzu))/(2*dz)

#     # this results in the following coupling coefficients:
#     return [sumeddyl+gravthermall, # down
#             sumeddyu-gravthermalu] # up
# end

# function fluxcoefs(z, dz, species, n_current)
#     # overload to generate the coefficients if they are not supplied
#     ntp = n_tot(n_current, z+dz)
#     nt0 = n_tot(n_current, z)
#     ntm = n_tot(n_current, z-dz)
#     Kp = Keddy(z+dz, ntp)
#     K0 = Keddy(z, nt0)
#     Km = Keddy(z-dz, ntm)
#     Tp = Temp(z+dz)
#     T0 = Temp(z)
#     Tm = Temp(z-dz)
#     Dp = Dcoef(Tp, ntp, species)
#     D0 = Dcoef(T0, nt0, species)
#     Dm = Dcoef(Tm, ntm, species)
#     Hsp = scaleH(z+dz, species)
#     Hs0 = scaleH(z, species)
#     Hsm = scaleH(z-dz, species)
#     H0p = scaleH(z+dz, Tp, n_current)
#     H00 = scaleH(z, T0, n_current)
#     H0m = scaleH(z-dz, Tm, n_current)

#     # return the coefficients
#     return fluxcoefs(z, dz,
#               [Km , K0, Kp],
#               [Dm , D0, Dp],
#               [Tm , T0, Tp],
#               [Hsm, Hs0, Hsp],
#               [H0m, H00, H0p],
#               species)
# end

# function lower_up(z, dz, species, n_current)
#     #= 
#     define transport coefficients for a given atmospheric layer for
#     transport from that layer to the one above. Variables ending in p refer to
#     the layer above ("plus"), in 0 refer to the current layer at altitude z,
#     and ending in m refer to the layer below ("minus") 

#     z: altitude
#     dz: altitude layer thickness ("resolution")
#     species: species symbol
#     n_current: current atmospheric state in a matrix by altitude and species.
#     =#
#     ntp = n_tot(n_current, z+dz)
#     nt0 = n_tot(n_current, z)
#     ntm = 1
#     Kp = Keddy(z+dz, ntp)
#     K0 = Keddy(z,nt0)
#     Km = 1
#     Tp = Temp(z+dz)
#     T0 = Temp(z)
#     Tm = 1
#     Dp = Dcoef(Tp, ntp, species)
#     D0 = Dcoef(T0, nt0, species)
#     Dm = 1
#     Hsp = scaleH(z+dz, species)
#     Hs0 = scaleH(z,species)
#     Hsm = 1
#     H0p = scaleH(z+dz, Tp, n_current)
#     H00 = scaleH(z,T0, n_current)
#     H0m = 1

#     # return the coefficients
#     return fluxcoefs(z, dz,
#               [Km , K0, Kp],
#               [Dm , D0, Dp],
#               [Tm , T0, Tp],
#               [Hsm, Hs0, Hsp],
#               [H0m, H00, H0p],
#               species)[2]
# end

# function upper_down(z, dz, species, n_current)
#     #= 
#     define transport coefficients for a given atmospheric layer for
#     transport from that layer to the one above. Variables ending in p refer to
#     the layer above ("plus"), in 0 refer to the current layer at altitude z,
#     and ending in m refer to the layer below ("minus") 

#     z: altitude
#     dz: altitude layer thickness ("resolution")
#     species: species symbol
#     n_current: current atmospheric state in a matrix by altitude and species.
#     =#
#     ntp = 1
#     nt0 = n_tot(n_current, z)
#     ntm = n_tot(n_current, z-dz)
#     Kp = 1
#     K0 = Keddy(z, nt0)
#     Km = Keddy(z-dz, ntm)
#     Tp = 1
#     T0 = Temp(z)
#     Tm = Temp(z-dz)
#     Dp = 1
#     D0 = Dcoef(T0, nt0, species)
#     Dm = Dcoef(Tm, ntm, species)
#     Hsp = 1
#     Hs0 = scaleH(z, species)
#     Hsm = scaleH(z-dz, species)
#     H0p = 1
#     H00 = scaleH(z, T0, n_current)
#     H0m = scaleH(z-dz, Tm, n_current)

#     # return the coefficients
#     return fluxcoefs(z, dz,
#               [Km , K0, Kp],
#               [Dm , D0, Dp],
#               [Tm , T0, Tp],
#               [Hsm, Hs0, Hsp],
#               [H0m, H00, H0p],
#               species)[1]
# end



# chemistry equations ==========================================================

# note: this whole section is basically witchcraft that Mike wrote and I've never
# really looked into it in detail. I just trust it works. --Eryn

function meanmass(n_current, z)
    #= find the mean molecular mass at a given altitude =#
    thisaltindex = n_alt_index[z]
    c = [n_current[sp][thisaltindex] for sp in specieslist]
    m = [speciesmolmasslist[sp] for sp in specieslist]
    return sum(c.*m)/sum(c)
end

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

# function reactionrates(n_current)
#     #=
#     Creates an array of size length(intaltgrid) x (number of reactions).
#     Populated with chemical reaction rates for each reaction based on species
#     populations.
#     =#
    
#     theserates = fill(convert(Float64, NaN),(length(intaltgrid),length(reactionnet)))
#     for ialt in 1:length(intaltgrid)
#         theserates[ialt,:] = reactionrates_local([[n_current[sp][ialt] for sp in specieslist];
#                                                 [n_current[J][ialt] for J in Jratelist];
#                                                 Temp(alt[ialt+1]);
#                                                 n_tot(n_current, alt[ialt+1])]...)
#     end
#     return theserates
# end

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

# plotting functions ===========================================================
function plotatm(n_current)
    clf()
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

    # selectspecies = [:O3, :H2O2, :HDO2, :DO2, :HO2, :HD, :H2, :H, :D, :H2O, 
    #                  :HDO, :O, :O2, :CO2]
    for sp in fullspecieslist
        plot(n_current[sp], alt[2:end-1]/1e5, color = speciescolor[sp],
             linewidth=2, label=sp, linestyle=speciesstyle[sp], zorder=1)
    end
    ylim(0, 250)
    xscale("log")
    xlim(1e-15, 1e18)
    xlabel(L"Species concentration (cm$^{-3}$)")
    ylabel("Altitude [km]")
    legend(bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
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




# boundary condition functions =================================================

function effusion_velocity(Texo::Float64, m::Float64, zmax)
    #=
    Returns effusion velocity for a species in cm/s
    Texo: temperature of the exobase (upper boundary) in K
    m: mass of one molecule of species in amu
    zmax: max altitude in cm
    =#
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+zmax))
    vth = sqrt(2*boltzmannK*Texo/(m*mH))  # this one is in m/s
    v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)  # this is in cm/s
    return v
end

function speciesbcs(species, surface_watersat, v_eff, oflux)
    H2Osat = surface_watersat["H2O"]
    HDOsat = surface_watersat["HDO"]

    speciesbclist=Dict(
                        :CO2=>["n" 2.1e17; "f" 0.],
                        :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                        :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                        :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                        :HDO=>["n" HDOsat[1]; "f" 0.],
                        :O=>["f" 0.; "f" oflux],
                        :H2=>["f" 0.; "v" v_eff["H2"]],
                        :HD=>["f" 0.; "v" v_eff["HD"]],
                        :H=>["f" 0.; "v" v_eff["H"]],
                        :D=>["f" 0.; "v" v_eff["D"]],
                      );
    get(speciesbclist, species, ["f" 0.; "f" 0.])
end

function boundaryconditions(species, dz, n_current, surf_watersat, v_eff, Of)
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

    bcs = speciesbcs(species, surf_watersat, v_eff, Of)
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

# TEMPERATURE ==================================================================

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

# WATER ========================================================================
# 1st term is a conversion factor to convert to (#/cm^3) bceause the 2nd
# term (from Washburn 1924) gives the value in mmHg
Psat(T::Float64) = ((133.3*1e-6)/(boltzmannK * T))*(10^(-2445.5646/T + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 - 6.757169))

# It doesn't matter to get the exact SVP of HDO because we never saturate. 
# I also tested it, and using the one I derived vs. the one for H2O makes no difference.
Psat_HDO(T::Float64) = ((133.3*1e-6)/(boltzmannK * T))*(10^(-2445.5646/T + 8.2312*log10(T) - 0.01677006*T + 1.20514e-5*T^2 - 6.757169))

# End module
end