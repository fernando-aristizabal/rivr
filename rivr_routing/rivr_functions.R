

# calculates friction slope from Manning's equation
friction_slope_mannings <- function(mannings_n, velocity, conversion_factor, hydraulic_radius) {

    friction_slope <- (  (mannings_n * velocity) / ( conversion_factor * (hydraulic_radius^(2/3)) )  )^2

    return(friction_slope)

}
