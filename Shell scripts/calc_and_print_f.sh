#!/bin/bash
# this file outputs the f value for converged scenarios. 
julia calc_f_converged.jl temp 190 110 200
julia calc_f_converged.jl temp 143 110 200
julia calc_f_converged.jl temp 238 110 200
julia calc_f_converged.jl temp 190 83 200
julia calc_f_converged.jl temp 190 138 200
julia calc_f_converged.jl temp 190 110 150
julia calc_f_converged.jl temp 190 110 250
julia calc_f_converged.jl water 2.1e-5
julia calc_f_converged.jl water 9.02e-3