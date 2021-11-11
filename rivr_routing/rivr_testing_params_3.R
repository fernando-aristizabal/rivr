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
max_times <- 25000
## initial flow condition (m^3/s) ##
initflows <- 4.020992
## conversion factor for Manning's equation ##
Cm <- 1.00
## gravitational acceleration (m/s^2) ##
g <- 9.81
## Set desired courant number for stability ##
# should be less than 0.06 for Dynamic wave
desired_courant_number <- 0.05
## output file template name ##
output_file_name <- '/data/rivr_outputs/outputs_20211108/rivr_output_params.h5'
## aggregate files ##
aggregate_files <- TRUE
## cores to use ##
cores <- 7
## boundary condition file ##
# set to NA to use native function # 
boundary_condition_file <- '/data/usgs_data/03602500_piney_river_vernon_tn_20210815_to_20210915.h5'
