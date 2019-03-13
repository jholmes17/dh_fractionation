#!/bin/bash
echo "Begin"
julia plot_Hflux.jl dh_3
julia plot_Hflux.jl dh_5.5
julia plot_Hflux.jl dh_8
julia plot_Hflux.jl dh_10
echo "Finished HD flux plots for DH experiments"
julia plot_Hflux.jl temp_192_83_199
julia plot_Hflux.jl temp_192_110_149
julia plot_Hflux.jl temp_192_110_199
julia plot_Hflux.jl temp_192_110_249
julia plot_Hflux.jl temp_192_138_199
julia plot_Hflux.jl temp_240_110_199
echo "Finished H flux plots for temp experiments"
julia plot_Hflux.jl water_1e-3
julia plot_Hflux.jl water_1e-4
julia plot_Hflux.jl water_1e-5
echo "Finished with everything"
