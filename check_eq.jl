################################################################################
# check_eq.jl
# TYPE: TEST
# WHICH: Equilibrium
# DESCRIPTION: Check if the atmosphere is in equilibrium by looking for
# Φ_O = 0.5Φ_{H+D}.
#
# Eryn Cangi
# 11 September 2018
# Last edited: 29 May 2020
# Currently tested for Julia: 1.4.1
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
main_results = results_dir * "VarWaterTemp/"
tradeoff_results = results_dir * "TradeoffPlots/Main Results/"
sim_rel_path = FNext * "/converged_" * FNext * ".h5"

# Set up the file path, automatically look for the simulation without external 
# specification of where it might be. 
# there may be duplicates in both folders but that's fine, it will only check once. 
if isfile(main_results * sim_rel_path)
    rf = main_results * sim_rel_path
elseif isfile(tradeoff_results * sim_rel_path)
    rf = tradeoff_results * sim_rel_path  
else
    throw("Invalid simulation parameters")
end



# Calculate the flux ratio =====================================================
Of = 1.2e8
Hf = get_flux(:H, rf, Of, temparr, therm_only=true)
HDf = Hf + get_flux(:D, rf, Of, temparr, therm_only=true)

# println("O flux: ", Of)
# println("H+D flux: ", HDf)
# println("H flux: ", Hf)

if round(HDf/Of) == 2
    println("Simulation: $(FNext)")
    println("Equilibrium: YES")
    println()
else
    println("Simulation: $(FNext)")
    println("Equilibrium: NO")
    println("Ratio: ", HDf/Of)
    println()
end


