#!/bin/bash
# run convergence for many temperature sims
# surface: 190-270, but lets go down to 150 per polar conditions from MCD
julia converge_new_file_250.jl temp 150 108 205 mean &
julia converge_new_file_250.jl temp 160 108 205 mean &
julia converge_new_file_250.jl temp 170 108 205 mean &
julia converge_new_file_250.jl temp 180 108 205 mean &
wait
julia converge_new_file_250.jl temp 190 108 205 mean &
julia converge_new_file_250.jl temp 200 108 205 mean &
julia converge_new_file_250.jl temp 210 108 205 mean &
julia converge_new_file_250.jl temp 216 108 205 mean & 
julia converge_new_file_250.jl temp 220 108 205 mean &
wait
julia converge_new_file_250.jl temp 230 108 205 mean &
julia converge_new_file_250.jl temp 240 108 205 mean &
julia converge_new_file_250.jl temp 250 108 205 mean &
julia converge_new_file_250.jl temp 260 108 205 mean &
wait
julia converge_new_file_250.jl temp 270 108 205 mean &

# tropopause: 75-160K
julia converge_new_file_250.jl temp 216 70 205 mean &
julia converge_new_file_250.jl temp 216 80 205 mean &
julia converge_new_file_250.jl temp 216 90 205 mean &
wait
julia converge_new_file_250.jl temp 216 100 205 mean &
julia converge_new_file_250.jl temp 216 110 205 mean &
julia converge_new_file_250.jl temp 216 120 205 mean &
julia converge_new_file_250.jl temp 216 130 205 mean &
wait
julia converge_new_file_250.jl temp 216 140 205 mean &
julia converge_new_file_250.jl temp 216 150 205 mean &
julia converge_new_file_250.jl temp 216 160 205 mean &


# exobase: 90-300K
# julia converge_new_file_250.jl temp 216 108 100 # this one makes no physical sense
julia converge_new_file_250.jl temp 216 108 125 mean &
wait
julia converge_new_file_250.jl temp 216 108 150 mean &
julia converge_new_file_250.jl temp 216 108 175 mean &
julia converge_new_file_250.jl temp 216 108 200 mean &
julia converge_new_file_250.jl temp 216 108 225 mean &
wait
julia converge_new_file_250.jl temp 216 108 250 mean &
julia converge_new_file_250.jl temp 216 108 275 mean &
julia converge_new_file_250.jl temp 216 108 300 mean &
julia converge_new_file_250.jl temp 216 108 325 mean &
wait
julia converge_new_file_250.jl temp 216 108 350 mean &