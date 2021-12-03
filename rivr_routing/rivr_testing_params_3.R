################################################################################
### Parameters to Use ###
################################################################################

## channel slope ( vertical m / horizontal m) ##
slopes <- c(0.001,0.0055,0.01)
## channel length (m) ##
extents <- 5000
## pde ##
pdes <- "Dynamic" #c("Kinematic","Dynamic")
## max spatial resolution (m) ##
max_dxs <- c(25,75,125)
## Manning's roughness ##
mannings <- c(0.1,0.15,0.20)
## channel bottom width (m) ## 
widths <- c(20,30,40)
## channel side slope (horizontal m / vertical m) ##
sideslopes <- c(0,1,2)
## time step (s) ##
dts <- 1
## max time step (s) ##
# this gets overwritten if maximum boundary file 
max_times <- 7200
## initial flow condition (m^3/s) ##
initflows <- 4 #4.020992
## conversion factor for Manning's equation ##
Cm <- 1.00
## gravitational acceleration (m/s^2) ##
g <- 9.81
## Set desired courant number for stability ##
# should be less than 0.06 for Dynamic wave
desired_courant_number <- 0.05
## output file template name ##
output_file_name <- '/data/rivr_outputs/3_outputs_testing/rivr_output_params.h5'
## aggregate files ##
aggregate_files <- TRUE
## cores to use ##
cores <- 15
## boundary condition file ##
# set to NA to use native function # 
#boundary_condition <- '/data/usgs_data/03602500_piney_river_vernon_tn_20210815_to_20210915.h5'
boundary_condition <- function(times,initflow) {

    num_of_times <- length(times)
    percentile_to_stop_wave <- 1/4
    index_to_stop_wave <- as.integer(percentile_to_stop_wave * num_of_times)
    time_to_stop_wave <- times[index_to_stop_wave]
    #times_wave <- times[1:index_to_stop_wave]
    #times_no_wave <- times[index_to_stop_wave+1,num_of_times]

    boundary_flows <- ifelse( 
                              times <= time_to_stop_wave,
                              initflow + (750/pi) * ( 1 - cos( pi * times / (percentile_to_stop_wave * index_to_stop_wave) ) ),
                              initflow
                            )
    #boundary_flows < - initflow + (750/pi) * ( 1 - cos( pi * times / num_of_times ) )
    
    return(boundary_flows)

}
