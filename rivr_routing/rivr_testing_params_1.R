################################################################################
### Parameters to Use ###
################################################################################

## channel slope ( vertical m / horizontal m) ##
slopes <- seq(0.001,0.02,0.003)
## channel length (m) ##
extents <- 1500
## Manning's roughness ##
mannings <- seq(0.04,0.20,0.02)
## channel bottom width (m) ## 
widths <- seq(20,50,5)
## channel side slope (horizontal m / vertical m) ##
sideslopes <- seq(0,1,0.2)
## time step (s) ##
dts <- 10
## max time step (s) ##
max_times <- 40000
## initial flow condition (m^3/s) ##
initflows <- seq(8,14,2)
## conversion factor for Manning's equation ##
Cm <- 1.00
## gravitational acceleration (m/s^2) ##
g <- 9.81
## Set desired courant number for stability ##
# should be less than 0.06 for Dynamic wave
desired_courant_number <- 0.04 
## output file template name ##
output_file_name <- '/data/outputs_20210608/rivr_output_params.h5'
## cores to use ##
cores <- 7
