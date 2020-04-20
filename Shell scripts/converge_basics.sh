#!/bin/bash

# temps - NEW
# Solar mean
# julia converge_new_file_250.jl temp 216 108 205 mean &
julia converge_new_file_250.jl temp 162 108 205 mean &
julia converge_new_file_250.jl temp 270 108 205 mean &
julia converge_new_file_250.jl temp 216 81 205 mean &
wait
julia converge_new_file_250.jl temp 216 135 205 mean &
julia converge_new_file_250.jl temp 216 108 158 mean &
julia converge_new_file_250.jl temp 216 108 264 mean &
julia converge_new_file_250.jl water 7.1e-7 mean &
wait
julia converge_new_file_250.jl water 1.33e-5 mean &
julia converge_new_file_250.jl water 1.41e-4 mean &
julia converge_new_file_250.jl water 3.72e-4 mean &
# julia converge_new_file_250.jl water 8.2e-4 mean &
wait
julia converge_new_file_250.jl water 2.89e-3 mean &

# Solar max
# julia converge_new_file_250.jl temp 216 108 205 smax &
# julia converge_new_file_250.jl temp 162 108 205 smax &
# julia converge_new_file_250.jl temp 270 108 205 smax &
# julia converge_new_file_250.jl temp 216 81 205 smax &
# wait
# julia converge_new_file_250.jl temp 216 135 205 smax &
# julia converge_new_file_250.jl temp 216 108 158 smax &
# julia converge_new_file_250.jl temp 216 108 264 smax &
# julia converge_new_file_250.jl water 1.48e-6 smax &
# wait
# julia converge_new_file_250.jl water 2.1e-5 smax &
# julia converge_new_file_250.jl water 8.1e-4 smax &
# julia converge_new_file_250.jl water 2.18e-3 smax &
# julia converge_new_file_250.jl water 4.46e-3 smax &
# wait
# julia converge_new_file_250.jl water 9.02e-3 smax &
# julia converge_new_file_250.jl water 1.358e-2 smax &

# Solar min
# julia converge_new_file_250.jl temp 216 108 205 smin &
# julia converge_new_file_250.jl temp 162 108 205 smin &
# julia converge_new_file_250.jl temp 270 108 205 smin &
# julia converge_new_file_250.jl temp 216 81 205 smin &
# wait
# julia converge_new_file_250.jl temp 216 135 205 smin &
# julia converge_new_file_250.jl temp 216 108 158 smin &
# julia converge_new_file_250.jl temp 216 108 264 smin &
# julia converge_new_file_250.jl water 1.48e-6 smin &
# wait
# julia converge_new_file_250.jl water 2.1e-5 smin &
# julia converge_new_file_250.jl water 8.1e-4 smin &
# julia converge_new_file_250.jl water 2.18e-3 smin &
# julia converge_new_file_250.jl water 4.46e-3 smin &
# wait
# julia converge_new_file_250.jl water 9.02e-3 smin &
# julia converge_new_file_250.jl water 1.358e-2 smin &

