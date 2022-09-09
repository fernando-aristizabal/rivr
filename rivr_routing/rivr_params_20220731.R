################################################################################
### Parameters to Use ###
################################################################################

## channel slope ( vertical m / horizontal m) ##
slopes <- c(0.005,0.0063,0.0076,0.0089,0.01)
## channel length (m) ##
extents <- 1000
## pde ##
pdes <- "Dynamic" #c("Kinematic","Dynamic")
## max spatial resolution (m) ##
dxs <- 10
## Manning's roughness ##
mannings <- seq(0.1,0.15,0.01)
## channel bottom width (m) ## 
widths <- seq(20,40,5)
## channel side slope (horizontal m / vertical m) ##
sideslopes <- 0
## time step (s) ##
dts <- 1
## max time step (s) ##
# this gets overwritten if maximum boundary file 
max_times <- 1000
## initial flow condition (m^3/s) ##
initflows <- c(4,6,8)
## monitoring node spacing ##
monitoring_node_step_sizes <- 1
## monitoring times spacing ##
monitoring_time_step_sizes <- 1
## conversion factor for Manning's equation ##
Cm <- 1.00
## gravitational acceleration (m/s^2) ##
g <- 9.81
## ioutput file template name ##
output_file_name <- '/data/rivr_outputs/outputs_20220731/rivr_output_params.h5'
## experiment id character length ##
experiment_id_character_length <- 5
## aggregate files ##
aggregate_files <- TRUE
## cores to use ##
cores <- 15
## boundary condition amplitude ##
amplitudes <- c(150,100,50)
## boundary condition cycles ##
numbers_of_cycles <- c(1,2)
## boundary condition wavelength as percentile of time ##
periods <- c(250,125)
## boundary condition function ##
boundary_condition <- function(times, initflow, amplitude, period, number_of_cycles) {

    # Fourier Series Article Wiki: https://en.wikipedia.org/wiki/Fourier_series

    #num_of_times <- length(times)
    #percentile_to_stop_wave <- 0.25
    #index_to_stop_wave <- as.integer(percentile_to_stop_wave * num_of_times)
    index_to_stop_wave <- which.min(abs(times - period))
    time_to_stop_wave <- times[index_to_stop_wave]
    times_wave <- times[1:index_to_stop_wave]
    
    # Koohafkan, Michael C., and Bassam A. Younis. "Open-channel computation with R." R Journal 7, no. 2 (2015): 249-262.
    boundary_flows <- ifelse( 
                             times <= time_to_stop_wave,
                             initflow + amplitude * ( 1 - cos( 2 * pi * number_of_cycles * times_wave / period ) ),
                             initflow
                            )
    
    return(boundary_flows)
    
}
