################################################################################
### Parameters to Use ###
################################################################################

## channel slope ( vertical m / horizontal m) ##
slopes <- c(0.005,0.008)
## channel length (m) ##
extents <- 10000
## pde ##
pdes <- "Dynamic" #c("Kinematic","Dynamic")
## max spatial resolution (m) ##
dxs <- 10
## Manning's roughness ##
mannings <- c(0.10,0.15)
## channel bottom width (m) ## 
widths <- c(40,60)
## channel side slope (horizontal m / vertical m) ##
sideslopes <- 0
## time step (s) ##
dts <- 1
## max time step (s) ##
# this gets overwritten if maximum boundary file 
max_times <- 1000
## initial flow condition (m^3/s) ##
initflows <- 4 #4.020992
## monitoring node spacing ##
monitoring_node_step_sizes <- 1
## monitoring times spacing ##
monitoring_time_step_sizes <- 1
## conversion factor for Manning's equation ##
Cm <- 1.00
## gravitational acceleration (m/s^2) ##
g <- 9.81
## output file template name ##
output_file_name <- '/data/rivr_outputs/outputs_20220310/rivr_output_params.h5'
## experiment id character length ##
experiment_id_character_length <- 4
## aggregate files ##
aggregate_files <- TRUE
## cores to use ##
cores <- 10
## boundary condition file ##
# set to NA to use native function # 
#boundary_condition <- '/data/usgs_data/03602500_piney_river_vernon_tn_20210815_to_20210915.h5'
boundary_condition <- function(times,initflow) {

    num_of_times <- length(times)
    percentile_to_stop_wave <- 0.25
    index_to_stop_wave <- as.integer(percentile_to_stop_wave * num_of_times)
    time_to_stop_wave <- times[index_to_stop_wave]
    times_wave <- times[1:index_to_stop_wave]

    boundary_flows <- ifelse( 
                              times <= time_to_stop_wave,
                              initflow + (500/pi) * ( 1 - cos( 2 * pi * times_wave / index_to_stop_wave ) ),
                              initflow
                            )
    #boundary_flows < - initflow + (750/pi) * ( 1 - cos( pi * times / num_of_times ) )
    
    return(boundary_flows)

}
