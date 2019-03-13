#!/bin/bash
echo "Running dh simulations"
julia -p 6 run_and_plot.jl dh 1
julia -p 6 run_and_plot.jl dh 3
julia -p 6 run_and_plot.jl dh 5.5
julia -p 6 run_and_plot.jl dh 8
julia -p 6 run_and_plot.jl dh 10
echo "finished dh simulations"

