#!/bin/bash
julia converge_new_file_250.jl water 1.48e-6 smean & # 0.1 pr μm
julia converge_new_file_250.jl water 2.1e-5 smean & # 1 pr μm
julia converge_new_file_250.jl water 8.1e-4 smean & # 10 pr μm
julia converge_new_file_250.jl water 2.18e-3 smean & # 25 pr μm
wait
julia converge_new_file_250.jl water 4.46e-3 smean & # 50 pr μm
julia converge_new_file_250.jl water 9.02e-3 smean & # 100 pr μm
julia converge_new_file_250.jl water 1.358e-2 smean & # 150 pr μm
