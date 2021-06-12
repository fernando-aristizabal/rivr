################################################################################
### Parameters to Use ###
################################################################################

## channel slope ( vertical m / horizontal m) ##
slopes <- c(0.009)
## channel length (m) ##
extents <- 1500
## Manning's roughness ##
mannings <- c(0.06,0.10,0.20)
## channel bottom width (m) ## 
widths <- c(20,30,50)
## channel side slope (horizontal m / vertical m) ##
sideslopes <- c(0,0.5,1)
## time step (s) ##
dts <- 2.5
## max time step (s) ##
max_times <- 30000
## initial flow condition (m^3/s) ##
initflows <- c(24,30,36)
## conversion factor for Manning's equation ##
Cm <- 1.00
## gravitational acceleration (m/s^2) ##
g <- 9.81
## Set desired courant number for stability ##
# should be less than 0.06 for Dynamic wave
desired_courant_number <- 0.05
## output file template name ##
output_file_name <- '/data/outputs_20210611/rivr_output_params.h5'
## cores to use ##
cores <- 7
