#!/bin/bash
echo "Running the simulations for the water variations"
julia -p 6 run_and_plot.jl water 1e-3
julia -p 6 run_and_plot.jl water 1e-4
julia -p 6 run_and_plot.jl water 1e-5
julia -p 6 run_and_plot.jl water 1e-6
echo "Finished water"

