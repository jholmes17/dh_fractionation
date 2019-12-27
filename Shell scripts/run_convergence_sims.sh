#!/bin/bash
echo "Beginning convergence script"
julia converge_new_file.jl temp 192 110 199 # new mean case
echo "Finished mean case"
echo "Starting cases where we lower some temperature"
julia converge_new_file.jl temp 144 110 199 # lower Tsurf
julia converge_new_file.jl temp 192 83 199 # lower Ttropo
julia converge_new_file.jl temp 192 110 149 # lower Texo
echo "Done with cases where we lower a temperature"
echo "Starting cases where we increase some temperature"
julia converge_new_file.jl temp 240 110 199 # higher Tsurf
julia converge_new_file.jl temp 192 138 199 # higher Ttropo
julia converge_new_file.jl temp 192 110 249 # higher Texo
echo "Finished temp cases"
# water cases...
echo "Starting water cases"
julia converge_new_file.jl water 1e-3
julia converge_new_file.jl water 1e-4
julia converge_new_file.jl water 1e-5
julia converge_new_file.jl water 1e-6
echo "Finished water cases"
echo "Starting DH cases"
julia converge_new_file.jl dh 1
julia converge_new_file.jl dh 3
julia converge_new_file.jl dh 5.5
julia converge_new_file.jl dh 8
julia converge_new_file.jl dh 10
echo "Finished DH cases"

echo "Completed all the convergences!"
