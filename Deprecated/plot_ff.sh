#!/bin/bash
#echo "Doing equator fractionation factor plots"
#julia plot_R_faster.jl eq-summer-morn
#julia plot_R_faster.jl eq-summer-afternoon
#julia plot_R_faster.jl eq-summer-night
#julia plot_R_faster.jl eq-winter-morn
#julia plot_R_faster.jl eq-winter-afternoon
#julia plot_R_faster.jl eq-winter-night
#echo "Finished with equator fractionation factor plots"

#echo "Doing North pole fractionation factor plots"
#julia plot_R_faster.jl n-pole-summer-morn
#julia plot_R_faster.jl n-pole-summer-afternoon
#julia plot_R_faster.jl n-pole-summer-night
#julia plot_R_faster.jl n-pole-winter-morn
#julia plot_R_faster.jl n-pole-winter-afternoon
#julia plot_R_faster.jl n-pole-winter-night
#echo "Finished with North pole fractionation factor plots"

echo "Doing South midlats fractionation factor plots"
julia plot_R_faster.jl s-midlat-summer-morn
julia plot_R_faster.jl s-midlat-summer-afternoon
julia plot_R_faster.jl s-midlat-summer-night
julia plot_R_faster.jl s-midlat-winter-morn
julia plot_R_faster.jl s-midlat-winter-afternoon
julia plot_R_faster.jl s-midlat-winter-night
echo "Finished with South midlats fractionation factor plots"

echo "Doing South pole fractionation factor plots"
julia plot_R_faster.jl s-pole-summer-morn
julia plot_R_faster.jl s-pole-summer-afternoon
julia plot_R_faster.jl s-pole-summer-night
#julia plot_R_faster.jl s-pole-winter-morn  # this simulation didn't work
julia plot_R_faster.jl s-pole-winter-afternoon
#julia plot_R_faster.jl s-pole-winter-night  # this simulation didn't work
echo "Finished with South pole fractionation factor plots"



