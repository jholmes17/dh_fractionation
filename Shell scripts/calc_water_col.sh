#!/bin/bash
echo "getting water col for temp simulations"
julia calc_water_col.jl temp 192 110 199 # new mean case
julia calc_water_col.jl temp 192 83 199 # lower Ttropo
julia calc_water_col.jl temp 192 110 149 # lower Texo
julia calc_water_col.jl temp 240 110 199 # higher Tsurf
julia calc_water_col.jl temp 192 138 199 # higher Ttropo
julia calc_water_col.jl temp 192 110 249 # higher Texo
echo "Finished temps"
echo "getting water col for water sims"
julia calc_water_col.jl water 1e-3
julia calc_water_col.jl water 1e-5
echo "Finished water"

