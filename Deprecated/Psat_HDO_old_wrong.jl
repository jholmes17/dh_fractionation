function Psat_HDO(T)
    #=
    Analytical expression for saturation vapor pressure of HDO in #/cm^3. The 
    conversion from N/m^2 to #/cm^3 is divide by N*m (that is, kT) and multiply 
    by the conversion to cm from m (1e-6). Justifications and sources are shown
    below.

    Input
        T: a single temperature in K
    Output:
        Pressure in #/cm^3
    =#
    function VP_h2o_over_liquid(T)
        #= 
        this subroutine is only used to find the VP of H2O at 217°C,
        the crossover point for HDO and H2O. The vapor pressure in hPa equation 
        comes from Jacobson eqn 2.61. 
        =# 
        T0C = 273.15 # temperature 0°C but in K
        VP_hPa = 6.112 * exp(6816*(1/T0C - 1/T) + 5.1309*log(T0C/T))
        return VP_hPa * (100/1) *(1e-6)/(boltzmannK*T)
    end
    
    # The following constants come from the equation for latent heat of 
    # sublimation (Jacobson eqn 2.56):
    # Ls = A + T(B + CT), valid for temps up to 273.15K. 
    A = 2.8345e6   # J/kg
    B = 340        # J/(kg K)
    C = 10.46      # J/(kg K^2)
    
    # This is the gas constant for HDO. from R = kB/m, where m is mass of species.
    # this equation found from Pierrehumbert's Principles of Planetary Climate.
    R = 434.8      # J/(kg K)
    
    # Reference pressure P0 at temperature T0. Here, T0 used is the crossover 
    # temperature for HDO and H2O, Jancso 1974, page 734. It's not perfect to 
    # use the crossover temp because it would be in the regime of VP over water 
    # rather than ice, but it's the best I can do. No other way for me to get a 
    # known VP for HDO, as far as I know.
    T0 = 273.15 + 217  # T0 in Kelvin
    P0 = VP_h2o_over_liquid(T0) * (boltzmannK*T/1) * (1/1e-6) * (1/100) # in hPa
    
    # The equation for VP in hPa below was found by integrating Jacobson eqn 
    # 2.63 using Jacobson eqn 2.56 for latent heat of sublimation. 
    VP_hPa = P0 * exp((1/R)*(A*(1/T0 - 1/T) + B*log(T0/T) + C*(T0-T))) # in hPa
    VP = VP_hPa * (100/1) *(1e-6)/(boltzmannK*T) # in #/cm^3.
    return VP
end