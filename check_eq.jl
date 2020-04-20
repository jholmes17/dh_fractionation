################################################################################
# check_eq.jl
# TYPE: TEST
# WHICH: Equilibrium
# DESCRIPTION: Check if the atmosphere is in equilibrium by looking for
# Φ_O = 0.5Φ_{H+D}.
#
# Eryn Cangi
# 11 September 2018
# Last edited: 3 January 2020
# Currently tested for Julia: 0.7
################################################################################
using HDF5
using Analysis

include("PARAMETERS.jl")

DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Experiment type and loading the file =========================================
# get arguments to use for accessing files 
argarray = Any[ARGS[i] for i in 1:1:length(ARGS)] 
if argarray[1] == "temp"
    FNext = "temp_$(argarray[2])_$(argarray[3])_$(argarray[4])"
    temparr = [parse(Float64, a) for a in argarray[2:end]]
else
    FNext = "$(argarray[1])_$(argarray[2])"
    temparr = [meanTs, meanTt, meanTe]
end

# load readfile and altitude array
lead = "/home/emc/GDrive-CU/Research/Results/"
rf = lead * FNext * "/converged_" * FNext * ".h5"

# Calculate the flux ratio =====================================================
Of = 1.2e8
Hf = get_flux(:H, "thermal", rf, Of, temparr)
HDf = Hf + get_flux(:D, "thermal", rf, Of, temparr)

println("O flux: ", Of)
println("H+D flux: ", HDf)
println("H flux: ", Hf)
println("Ratio: ", HDf/Of)
println()
