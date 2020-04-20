#!/bin/bash
# run convergence for many temperature sims
# surface: 190-270, but lets go down to 150 per polar conditions from MCD
julia converge_new_file_250.jl temp 150 130 205 mean &
julia converge_new_file_250.jl temp 160 130 205 mean &
julia converge_new_file_250.jl temp 170 130 205 mean &
julia converge_new_file_250.jl temp 180 130 205 mean &
wait
julia converge_new_file_250.jl temp 190 130 205 mean &
julia converge_new_file_250.jl temp 200 130 205 mean &
julia converge_new_file_250.jl temp 210 130 205 mean &
julia converge_new_file_250.jl temp 220 130 205 mean &
wait
julia converge_new_file_250.jl temp 230 130 205 mean &
julia converge_new_file_250.jl temp 240 130 205 mean &
julia converge_new_file_250.jl temp 250 130 205 mean &
julia converge_new_file_250.jl temp 260 130 205 mean &
wait
julia converge_new_file_250.jl temp 270 130 205 mean &

# tropopause: 75-160K
julia converge_new_file_250.jl temp 216 100 205 mean &
julia converge_new_file_250.jl temp 216 110 205 mean &
julia converge_new_file_250.jl temp 216 120 205 mean &
wait
julia converge_new_file_250.jl temp 216 130 205 mean &
julia converge_new_file_250.jl temp 216 140 205 mean &
julia converge_new_file_250.jl temp 216 150 205 mean &
julia converge_new_file_250.jl temp 216 160 205 mean &
wait

# exobase: 90-300K
julia converge_new_file_250.jl temp 216 130 150 mean &
julia converge_new_file_250.jl temp 216 130 175 mean &
julia converge_new_file_250.jl temp 216 130 200 mean &
julia converge_new_file_250.jl temp 216 130 225 mean &
wait
julia converge_new_file_250.jl temp 216 130 250 mean &
julia converge_new_file_250.jl temp 216 130 275 mean &
julia converge_new_file_250.jl temp 216 130 300 mean &
julia converge_new_file_250.jl temp 216 130 325 mean &
wait
julia converge_new_file_250.jl temp 216 130 350 mean &