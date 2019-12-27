#!/bin/bash
echo "Begin"
julia plot_FF.jl dh_1
julia plot_FF.jl dh_3
julia plot_FF.jl dh_5.5
julia plot_FF.jl dh_8
julia plot_FF.jl dh_10
echo "Finished plots for DH experiments"
julia plot_FF.jl temp_192_110_199 # mean
julia plot_FF.jl temp_192_110_149 # vary exobase
julia plot_FF.jl temp_192_110_249 
julia plot_FF.jl temp_192_83_199 # vary tropopause
julia plot_FF.jl temp_192_138_199
julia plot_FF.jl temp_240_110_199 # vary surface
echo "Finished plots for temp experiments"
julia plot_FF.jl water_1e-3
julia plot_FF.jl water_1e-4
julia plot_FF.jl water_1e-5
echo "Finished with everything"
