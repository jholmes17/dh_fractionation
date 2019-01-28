## run_and_plot.jl --- routines to run the coupled photochemistry model
## to produce model output useful for a 2016 photochemistry paper

arg_from_rnp = Any[ARGS[i] for i in 1:1:length(ARGS)] # get the args from command line
@eval @everywhere arg_from_rnp=$arg_from_rnp          # make accessible to all scripts

# following code is specifically for temp/water/dh variation experiments
# extfn = extension to a filename - where a specific experiment is
if arg_from_rnp[1] == "temp"
    extfn = "temp_$(arg_from_rnp[2])_$(arg_from_rnp[3])_$(arg_from_rnp[4])"
elseif arg_from_rnp[1] == "water"
    extfn = "water_$(arg_from_rnp[2])"
elseif arg_from_rnp[1] == "dh"
    extfn = "dh_$(arg_from_rnp[2])"
end

@eval @everywhere extfn=$extfn
@everywhere include("setup_photochemistry.jl") # NOTE: change this as needed
@everywhere @time update!(n_current,0.)

# Run the simulation with logarithmic time steps
@everywhere timepts = logspace(log10(1),log10(1e7*3.14e7),1000)
@everywhere timediff = timepts[2:end]-timepts[1:end-1]
@everywhere append!(timediff,3.3e6*3.14e7*ones(Float64,300))

# identify the converged test case file NOTE: fix all this as needed
scriptdir = "/data/GoogleDrive/Phys/LASP/chaffincode-working/" # the main script directory
experimentdir = "/data/VaryTW_Ana/"  # used for when experiments stored on desktop
# readfile = scriptdir*extfn*"/converged_standardwater_D_"*extfn*".h5"  # if experiments stored in google drive
readfile = experimentdir*extfn*"/converged_standardwater_D_"*extfn*".h5"
println("ALERT: Using file: ", readfile)

# water ppm of 20, 40, 60, 80 each tested at altitudes 20, 40, 60, 80, 120
# set up the pairings and some filenames
waterppmvec = [20 40 60 80]
hdoppmvec = waterppmvec * DH
wateraltvec = [20 40 60 80 100 120]
parmsvec = [[a,b] for a in waterppmvec, b in wateraltvec]
parmsvec = reshape(parmsvec,length(parmsvec))
filenamevec = [string(experimentdir*extfn*"/ppm_", a[1], "_alt_", a[2], ".h5")
               for a in parmsvec]

@everywhere function runwaterprofile(n_current, ppmadd, peakalt, dtlist, filename)\
    #=
    n_current: matrix of concentrations by species and altitude, with water
               profile modified to include an enhanced water concentration of
               a certain ppm at a certain altitudes
    ppmadd:    number representing the enhanced ppm of water to use
    peakalt:   number representing the altitude at which the water parcel is
               added
    dtlist:    list of times at which to recalculate the matrix of concentrations
    fliename:  filename specifying ppm and altitude for results to be written to

    Updates the n_current to include the enhanced water ppm at a certain
    altitude.
    =#

    println("Now working on ppm $(ppmadd) and alt $(peakalt)")
    n_internal = deepcopy(n_current)

    # modify the water profile stored in n_internal with "gaussian" packet
    waterppm = 1e-6*map(x->ppmadd.*exp(-((x-peakalt)/12.5)^2),alt[2:end-1]/1e5)+H2Oinitfrac
    waterprofile = waterppm.*map(z->n_tot(n_internal,z),alt[2:end-1])
    n_internal[:H2O] = waterprofile

    # modify the HDO profile stored in n_internal with "gaussian" packet
    ppmadd_HDO = ppmadd*DH
    hdoppm = 1e-6*map(x->ppmadd_HDO.*exp(-((x-peakalt)/12.5)^2),alt[2:end-1]/1e5)+HDOinitfrac
    hdoprofile=hdoppm.*map(z->n_tot(n_internal,z),alt[2:end-1])
    n_internal[:HDO] = hdoprofile

    println("About to run the $(ppmadd) ppm and $(peakalt) alt profile")
    result = runprofile(n_internal, dtlist, filename)
    if result==nothing
        println("skipped $(ppmadd) ppm / alt $(peakalt)")
    else
        println("Successfully finished $(ppmadd) ppm and alt $(peakalt)")
    end
    println()
    return result
end

@everywhere function runprofile(n_current, dtlist, filename)
    #=
    n_current: matrix of concentrations by species and altitude, with water
               profile modified to include an enhanced water concentration of
               a certain ppm at a certain altitudes
    dtlist:    list of times at which to find the new concentration matrix
    filename:  filename specifying ppm concentration and altitude of the parcel
               to which the results will be written

    This function actually does the work of running the simulation and returning
    the resulting matrix of concentrations by species and altitude.
    =#
    n_internal = deepcopy(n_current)
    elapsed_time = 0.0

    # create a matrix to contain the data in n_internal, which is a Dict
    n_internal_mat = Array{Float64}(length(alt)-2,length(collect(keys(n_internal))));
    for ispecies in 1:length(collect(keys(n_internal)))
        for ialt in 1:length(alt)-2
            n_internal_mat[ialt,ispecies] = n_internal[collect(keys(n_internal))[ispecies]][ialt]
        end
    end

    if isfile(filename)==true
        return nothing
    elseif isfile(filename)==false
        h5write(filename,"n_current/init",n_internal_mat)
        h5write(filename,"n_current/alt",alt)
        h5write(filename,"n_current/species",map(string,collect(keys(n_internal))))
        h5write(filename,"n_current/timelist",cumsum(dtlist))

        # TIME LOOP - simulate over time
        thisi=0
        for dt in dtlist
            #println(filename*": iteration = "* string(thisi+=1)*" "*Libc.strftime(time()))
            thisi += 1
            elapsed_time+=dt
            ##n_old=deepcopy(n_internal)  #TODO: what's this line? remove it?

            # where the action happens - update n_internal for new timesteps
            update!(n_internal,dt)

            # save the concentrations to history and write n_current into
            # n_current_mat
            n_internal_mat = Array{Float64}(length(alt)-2,length(collect(keys(n_internal))));
            for ispecies in 1:length(collect(keys(n_internal)))
                for ialt in 1:length(alt)-2
                    n_internal_mat[ialt,ispecies] = n_internal[collect(keys(n_internal))[ispecies]][ialt]
                end
            end

            # write n_internal_mat to file
            h5write(filename, string("n_current/iter_",thisi), n_internal_mat)
        end
        return n_internal
    end
end

@everywhere function read_ncurrent_from_file(readfile,tag)
    thisalt = h5read(readfile,"n_current/alt")
    if thisalt != alt
        throw("altitudes in file do not match altitudes in memory!")
    end
    n_current_tag_list = map(Symbol,h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,tag);
    n_current = Dict{Symbol,Array{Float64,1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]]=reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

@everywhere function get_H_fluxes(readfile)
    #=
    Produces an array of H fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
    =#
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)
    Hfluxes = fill(0.,timelength)
    n_current = read_ncurrent_from_file(readfile,string("n_current/init"))

    # Calculate the total H flux: sum of (H population @ 200 km) * (H flux rate)
    # and (H_2 population @ 200 km) * (H_2 flux rate )
    Hfluxes[1] = (n_current[:H][end]*speciesbcs(:H)[2,2]
                  + 2*n_current[:H2][end]*speciesbcs(:H2)[2,2]
                  + n_current[:HD][end]*speciesbcs(:HD)[2,2])
    for i in 1:(timelength-1)
        n_current=read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        Hfluxes[i+1]=(n_current[:H][end]*speciesbcs(:H)[2,2]
                      + 2*n_current[:H2][end]*speciesbcs(:H2)[2,2]
                      + n_current[:HD][end]*speciesbcs(:HD)[2,2])
    end
    return Hfluxes
end

@everywhere function get_D_fluxes(readfile)
    #=
    Produces an array of D fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
    =#
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)
    Dfluxes = fill(0.,timelength)
    n_current = read_ncurrent_from_file(readfile,string("n_current/init"))

    # Calculate the total H flux: sum of (H population @ 200 km) * (H outward flux)
    # and (H_2 population @ 200 km) * (H_2 outward flux)
    Dfluxes[1] = (n_current[:D][end]*speciesbcs(:D)[2,2]
                  +n_current[:HD][end]*speciesbcs(:HD)[2,2])
    for i in 1:(timelength-1)
        n_current=read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        Dfluxes[i+1]=(n_current[:D][end]*speciesbcs(:D)[2,2]
                      +n_current[:HD][end]*speciesbcs(:HD)[2,2])
    end
    return Dfluxes
end

@everywhere function get_H_and_D_fluxes(readfile)
    #=
    Produces an array of H and D fluxes at the top of the atmosphere (200 km) due to
    introduction of water parcels of various ppm and various altitudes.
    =#
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)
    HDfluxes = fill(0.,timelength)
    n_current = read_ncurrent_from_file(readfile,string("n_current/init"))

    # Calculate the total H+D flux: sum of (H population @ 200 km) * (H flux rate)
    # and (H_2 population @ 200 km) * (H_2 flux rate), and deuterated species
    HDfluxes[1] = (n_current[:H][end]*speciesbcs(:H)[2,2]       # H
                  + 2*n_current[:H2][end]*speciesbcs(:H2)[2,2]  # H
                  + n_current[:HD][end]*speciesbcs(:HD)[2,2]    # H
                  + n_current[:D][end]*speciesbcs(:D)[2,2]      # D
                  + n_current[:HD][end]*speciesbcs(:HD)[2,2]) # D
    for i in 1:(timelength-1)
        n_current=read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        HDfluxes[i+1] = (n_current[:H][end]*speciesbcs(:H)[2,2]
                      + 2*n_current[:H2][end]*speciesbcs(:H2)[2,2]
                      + n_current[:HD][end]*speciesbcs(:HD)[2,2]
                      + n_current[:D][end]*speciesbcs(:D)[2,2]
                      + n_current[:HD][end]*speciesbcs(:HD)[2,2])
    end
    return HDfluxes
end

@everywhere function get_rates_and_fluxes(readfile)
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)
    reactionrateshist = fill(convert(Float64,NaN),timelength,length(intaltgrid),length(reactionnet))
    fluxhist = fill(convert(Float64,NaN),timelength,length(intaltgrid),length(specieslist))
    n_current = read_ncurrent_from_file(readfile,string("n_current/init"))
    reactionrateshist[1,:,:] = reactionrates(n_current)
    fluxhist[1,:,:] = fluxes(n_current,dz)
    println("doing the loop")
    for i in 1:(timelength-1)
        n_current = read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        reactionrateshist[i+1,:,:] = reactionrates(n_current)
        fluxhist[i+1,:,:] = fluxes(n_current,dz)
    end
    println("done with loop")
    (reactionrateshist,fluxhist)
end

@everywhere function get_all_rates_and_fluxes(readfile)
    println("going into get_rates_and_fluxes")
    (reactionrateshist,fluxhist)=get_rates_and_fluxes(readfile)

    # doing this way makes it work if files already exist
    h5open(readfile, isfile(readfile) ? "r+" : "w") do file
       write(file,"fluxes/flux_history",fluxhist)
       write(file,"rates/reaction_rates_history",reactionrateshist)
    end
    return
end

# Functions to retrieve water profiles
function get_water_ppm(filename)
    n_file = read_ncurrent_from_file(filename,"n_current/init")
    waterppmvec = 1e6*n_file[:H2O]./map(z->n_tot(n_file,z),alt[2:end-1])
    return waterppmvec
end

function get_hdo_ppm(filename)
    n_file = read_ncurrent_from_file(filename,"n_current/init")
    hdoppmvec = 1e6*n_file[:HDO]./map(z->n_tot(n_file,z),alt[2:end-1])
    return hdoppmvec
end

pmap(x->println(string("parmsvec[i][1]=",x[1],", parmsvec[i][2]=",x[2],",
     filename=",x[3])),[[p,f;] for (p,f) in zip(parmsvec,filenamevec)])

# This runs the simulation for a year and returns
oneyeartimepts = logspace(log10(1),log10(3.14e7),1000)
oneyeartimediff = oneyeartimepts[2:end]-oneyeartimepts[1:end-1]
n_converged = get_ncurrent(readfile)

# Add water / run for a year / remove water / run for a year ===================
println("running sim for one year")
oneyearfn = experimentdir*extfn*"/one_year_response_to_80ppm_at_60km.h5"
n_oneyear = runwaterprofile(n_converged, 80, 60, oneyeartimediff, oneyearfn)

println("now removing the water")
returnfn = experimentdir*extfn*"/one_year_response_to_80ppm_at_60km_return.h5"
n_return = runwaterprofile(n_oneyear, 0., 60, oneyeartimediff, returnfn)

# Add water and run for just over a year, no removal ===========================
# This runs the simulation for all added ppms and altitudes
println("Now doing water profiles")
pmap(x->runwaterprofile(n_current, x[1], x[2], timediff, x[3]),
     [[p,f;] for (p,f) in zip(parmsvec,filenamevec)])
println("Finished with water profiles")
pmap(get_all_rates_and_fluxes,filenamevec)

# Calculate H, D, H+D flux at exobase due to each experiment ===================
println("Doing H fluxes")
Hfluxes = pmap(get_H_fluxes,filenamevec)
lhfl = length(Hfluxes[1,1])
writeHfluxes = fill(0.0,(length(waterppmvec),length(wateraltvec),lhfl+2))
for lp in 1:length(parmsvec)
    ippm = lp%length(waterppmvec)+1
    ialt = floor(Int,(lp-1)/length(waterppmvec))+1
    writeHfluxes[ippm,ialt,:] = [parmsvec[lp],Hfluxes[lp];]
end

println("Doing D fluxes")
Dfluxes = pmap(get_D_fluxes,filenamevec)
lhfl = length(Dfluxes[1,1])
writeDfluxes = fill(0.0,(length(waterppmvec),length(wateraltvec),lhfl+2))
for lp in 1:length(parmsvec)
    ippm = lp%length(waterppmvec)+1
    ialt = floor(Int,(lp-1)/length(waterppmvec))+1
    writeDfluxes[ippm,ialt,:] = [parmsvec[lp],Dfluxes[lp];]
end

println("Doing H+D fluxes")
HDfluxes = pmap(get_H_and_D_fluxes,filenamevec)
lhfl = length(HDfluxes[1,1])
writeHDfluxes = fill(0.0,(length(hdoppmvec),length(wateraltvec),lhfl+2))
for lp in 1:length(parmsvec)
    ippm = lp%length(hdoppmvec)+1
    ialt = floor(Int,(lp-1)/length(hdoppmvec))+1
    writeHDfluxes[ippm,ialt,:] = [parmsvec[lp],HDfluxes[lp];]
end

# Calculate the water profiles =================================================
waterprofs = map(get_water_ppm,filenamevec)
writewaterprof = fill(0.0,(length(waterppmvec),length(wateraltvec),length(alt)-2+2))
hdoprofs = map(get_hdo_ppm,filenamevec)
writehdoprof = fill(0.0,(length(hdoppmvec),length(wateraltvec),length(alt)-2+2))
for lp in 1:length(parmsvec)
    ippm = lp%length(waterppmvec)+1
    ialt = floor(Int,(lp-1)/length(waterppmvec))+1
    writewaterprof[ippm,ialt,:] = [parmsvec[lp],waterprofs[lp];]
    writehdoprof[ippm,ialt,:] = [parmsvec[lp],hdoprofs[lp];]
end


# WRITE OUT H, D, H+D ESCAPE FLUX HISTORY ======================================
println("Writing H esc file")
hfile = experimentdir*extfn*"/H_esc_flux_history.h5"
h5open(hfile, isfile(hfile) ? "r+" : "w") do file
   write(file,"fluxes/fluxvals",writeHfluxes)
   write(file,"fluxes/times",h5read(experimentdir*extfn*"/ppm_20_alt_20.h5","n_current/timelist"))
   write(file,"waterprofs/ppm",writewaterprof)
   write(file,"waterprofs/alt",alt[2:end-1])
end

println("Writing D esc file")
hfile = experimentdir*extfn*"/D_esc_flux_history.h5"
h5open(hfile, isfile(hfile) ? "r+" : "w") do file
   write(file,"fluxes/fluxvals",writeDfluxes)
   write(file,"fluxes/times",h5read(experimentdir*extfn*"/ppm_20_alt_20.h5","n_current/timelist"))
   write(file,"hdoprofs/ppm",writehdoprof)
   write(file,"waterprofs/alt",alt[2:end-1])
end

println("Writing H+D esc file")
hdfile = experimentdir*extfn*"/H_and_D_esc_flux_history.h5"
h5open(hdfile, isfile(hdfile) ? "r+" : "w") do file
   write(file,"fluxes/fluxvals",writeHDfluxes)
   write(file,"fluxes/times",h5read(experimentdir*extfn*"/ppm_20_alt_20.h5","n_current/timelist"))
   write(file,"waterprofs/ppm",writewaterprof)
   write(file, "hdoprofs/ppm", writehdoprof)
   write(file,"waterprofs/alt",alt[2:end-1])
end
