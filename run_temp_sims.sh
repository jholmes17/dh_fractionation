#!/bin/bash
echo "Running temp simulations"
julia -p 6 run_and_plot.jl temp 192 110 199 # new mean case
julia -p 6 run_and_plot.jl temp 144 110 199 # lower Tsurf
julia -p 6 run_and_plot.jl temp 192 83 199 # lower Ttropo
julia -p 6 run_and_plot.jl temp 192 110 149 # lower Texo
julia -p 6 run_and_plot.jl temp 240 110 199 # higher Tsurf
julia -p 6 run_and_plot.jl temp 192 138 199 # higher Ttropo
julia -p 6 run_and_plot.jl temp 192 110 249 # higher Texo
echo "Finished temps"
