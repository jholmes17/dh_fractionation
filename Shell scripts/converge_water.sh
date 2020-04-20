#!/bin/bash
# NOTICE: THESE NUMBERS ARE OLD - DO NOT USE THIS FILE
# JUST USE THE WATER EXPERIMENT RESULTS FROM "CONVERGE_BASICS"
julia converge_new_file_250.jl water 2.1e-5 mean & # 1 pr μm
julia converge_new_file_250.jl water 8.1e-4 mean & # 10 pr μm
julia converge_new_file_250.jl water 2.18e-3 mean & # 25 pr μm
julia converge_new_file_250.jl water 4.46e-3 mean & # 50 pr μm
wait
julia converge_new_file_250.jl water 9.02e-3 mean & # 100 pr μm

