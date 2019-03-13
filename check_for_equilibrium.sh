#!/bin/bash
echo "Checking equilibrium"
julia check_eq.jl temp_192_83_199 temp 192 83 199
julia check_eq.jl temp_192_110_149 temp 192 110 149
julia check_eq.jl temp_192_110_199 temp 192 110 199
julia check_eq.jl temp_192_138_199 temp 192 138 199
julia check_eq.jl temp_192_110_249 temp 192 110 249
julia check_eq.jl temp_240_110_199 temp 240 110 199


echo "Checking equilibrium, water"
julia check_eq.jl water_1e-3 water
julia check_eq.jl water_1e-4 water
julia check_eq.jl water_1e-5 water
echo "FINISHED"


