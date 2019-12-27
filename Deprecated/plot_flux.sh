#!/bin/bash
#echo "Starting North pole flux plots"
#julia plot_Hflux.jl n-pole-summer-morn
#julia plot_Hflux.jl n-pole-summer-afternoon
#julia plot_Hflux.jl n-pole-summer-night
#julia plot_Hflux.jl n-pole-winter-morn
#julia plot_Hflux.jl n-pole-winter-afternoon
#julia plot_Hflux.jl n-pole-winter-night
#echo "Finished with North pole flux plots"

echo "Starting North midlats flux plots"
julia plot_Hflux.jl n-midlat-summer-morn
julia plot_Hflux.jl n-midlat-summer-afternoon
julia plot_Hflux.jl n-midlat-summer-night
julia plot_Hflux.jl n-midlat-winter-morn
julia plot_Hflux.jl n-midlat-winter-afternoon
julia plot_Hflux.jl n-midlat-winter-night
echo "Finished with North midlats flux plots"

#echo "Starting South midlats flux plots"
#julia plot_Hflux.jl s-midlat-summer-morn
#julia plot_Hflux.jl s-midlat-summer-afternoon
#julia plot_Hflux.jl s-midlat-summer-night
#julia plot_Hflux.jl s-midlat-winter-morn
#julia plot_Hflux.jl s-midlat-winter-afternoon
#julia plot_Hflux.jl s-midlat-winter-night
#echo "Finished with South midlats flux plots"

#echo "Starting South pole flux plots"
#julia plot_Hflux.jl s-pole-summer-morn
#julia plot_Hflux.jl s-pole-summer-afternoon
#julia plot_Hflux.jl s-pole-summer-night
#julia plot_Hflux.jl s-pole-winter-morn  # this simulation didn't work
#julia plot_Hflux.jl s-pole-winter-afternoon
#julia plot_Hflux.jl s-pole-winter-night  # this simulation didn't work
#echo "Finished with South pole flux plots"

#!/bin/bash
#echo "Starting equator flux plots"
#julia plot_Hflux.jl eq-summer-morn
#julia plot_Hflux.jl eq-summer-afternoon
#julia plot_Hflux.jl eq-summer-night
#julia plot_Hflux.jl eq-winter-morn
#julia plot_Hflux.jl eq-winter-afternoon
#julia plot_Hflux.jl eq-winter-night
#echo "Finished with equator flux plots"




