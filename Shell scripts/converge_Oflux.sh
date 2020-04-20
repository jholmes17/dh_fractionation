#!/bin/bash
# run convergence for many O flux values
# julia converge_new_file_250.jl Oflux 3e7 mean &
# julia converge_new_file_250.jl Oflux 4e7 mean &
# julia converge_new_file_250.jl Oflux 5e7 mean &
# julia converge_new_file_250.jl Oflux 6e7 mean &
# wait
# julia converge_new_file_250.jl Oflux 7e7 mean &
julia converge_new_file_250.jl Oflux 8e7 mean &
julia converge_new_file_250.jl Oflux 9e7 mean &
julia converge_new_file_250.jl Oflux 1.0e8 mean &
wait
julia converge_new_file_250.jl Oflux 1.1e8 mean &
julia converge_new_file_250.jl Oflux 1.2e8 mean &
julia converge_new_file_250.jl Oflux 1.3e8 mean &
julia converge_new_file_250.jl Oflux 1.4e8 mean &
wait
julia converge_new_file_250.jl Oflux 1.5e8 mean &
julia converge_new_file_250.jl Oflux 1.6e8 mean &