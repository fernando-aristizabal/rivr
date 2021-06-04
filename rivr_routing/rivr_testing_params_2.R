################################################################################
### Parameters to Use ###
################################################################################

## channel slope ( vertical m / horizontal m) ##
slopes <- 0.001
## channel length (m) ##
extents <- 1500
## Manning's roughness ##
mannings <- c(0.04,0.06,0.08)
## channel bottom width (m) ## 
widths <- 20
## channel side slope (horizontal m / vertical m) ##
sideslopes <- 0
## time step (s) ##
dts <- 10
## max time step (s) ##
max_times <- 40000
## initial flow condition (m^3/s) ##
initflows <- 6
## conversion factor for Manning's equation ##
Cm <- 1.00
## gravitational acceleration (m/s^2) ##
g <- 9.81
## Set desired courant number for stability ##
# should be less than 0.06 for Dynamic wave
desired_courant_number <- 0.04 
## output file template name ##
output_file_name <- '/data/outputs_20210526/rivr_output_params.h5'
## cores to use ##
cores <- 3
