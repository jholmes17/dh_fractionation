#!/bin/bash
echo "Checking equilibrium"
julia check_eq.jl temp 192 83 199
julia check_eq.jl temp 192 110 149
julia check_eq.jl temp 192 110 199
julia check_eq.jl temp 192 138 199
julia check_eq.jl temp 192 110 249
julia check_eq.jl temp 240 110 199


echo "Checking equilibrium, water"
julia check_eq.jl water 1e-3
julia check_eq.jl water 1e-4
julia check_eq.jl water 1e-5
echo "FINISHED"


