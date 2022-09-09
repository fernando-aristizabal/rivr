library(stringr)


# calculates friction slope from Manning's equation
friction_slope_mannings <- function(Q, C) {

    #friction_slope <- (  (mannings_n * velocity) / ( conversion_factor * (hydraulic_radius^(2/3)) )  )^2
    friction_slope <- (Q/C) ^ 2.0

    return(friction_slope)

}


# converts experiment numbers to ids'
convert_experiment_number_to_ids <- function(experiment_number, id_length) {
    
    experiment_id <- stringr::str_pad(experiment_number,id_length,"0",side='left')

    return(experiment_id)
}
    
    
    
    # trying to find optimal spatial resolution, dx, to optimize courant number
    #max_velocity <- -10e6
    #for (b in boundary) {
    #  norm_dep <- rivr::normal_depth(slope, manning, b, 1, Cm, width, sideslope)
    #  channel_geometry <- rivr::channel_geom(norm_dep,width,sideslope)
    #  velocity <- b / channel_geometry[['A']]
    #  if (velocity > max_velocity) {
    #      max_velocity <- velocity
    #  }
    #}
    #dx <- max_velocity * dt / desired_courant_number

    # sets dx to maximum dx value desired if not already
    #if (dx > max_dx) {
    #      dx <- max_dx
    #}
