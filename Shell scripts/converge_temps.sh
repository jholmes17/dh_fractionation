#!/bin/bash
# run convergence for many temperature sims
# surface: 190-270, but lets go down to 150 per polar conditions from MCD
julia converge_new_file_250.jl temp 150 108 205 smean &
julia converge_new_file_250.jl temp 160 108 205 smean &
julia converge_new_file_250.jl temp 170 108 205 smean &
julia converge_new_file_250.jl temp 180 108 205 smean &
wait
julia converge_new_file_250.jl temp 190 108 205 smean &
julia converge_new_file_250.jl temp 200 108 205 smean &
julia converge_new_file_250.jl temp 210 108 205 smean &
julia converge_new_file_250.jl temp 220 108 205 smean &
wait
julia converge_new_file_250.jl temp 230 108 205 smean &
julia converge_new_file_250.jl temp 240 108 205 smean &
julia converge_new_file_250.jl temp 250 108 205 smean &
julia converge_new_file_250.jl temp 260 108 205 smean &
wait
julia converge_new_file_250.jl temp 270 108 205 smean &

# tropopause: 75-160K
julia converge_new_file_250.jl temp 216 70 205 smean &
julia converge_new_file_250.jl temp 216 80 205 smean &
julia converge_new_file_250.jl temp 216 90 205 smean &
wait
julia converge_new_file_250.jl temp 216 100 205 smean &
julia converge_new_file_250.jl temp 216 110 205 smean &
julia converge_new_file_250.jl temp 216 120 205 smean &
julia converge_new_file_250.jl temp 216 130 205 smean &
wait
julia converge_new_file_250.jl temp 216 140 205 smean &
julia converge_new_file_250.jl temp 216 150 205 smean &
julia converge_new_file_250.jl temp 216 160 205 smean &


# exobase: 90-300K
# julia converge_new_file_250.jl temp 216 108 100 # this one makes no physical sense
julia converge_new_file_250.jl temp 216 108 125 smean &
wait
julia converge_new_file_250.jl temp 216 108 150 smean &
julia converge_new_file_250.jl temp 216 108 175 smean &
julia converge_new_file_250.jl temp 216 108 200 smean &
julia converge_new_file_250.jl temp 216 108 225 smean &
wait
julia converge_new_file_250.jl temp 216 108 250 smean &
julia converge_new_file_250.jl temp 216 108 275 smean &
julia converge_new_file_250.jl temp 216 108 300 smean &
julia converge_new_file_250.jl temp 216 108 325 smean &
wait
julia converge_new_file_250.jl temp 216 108 350 smean &