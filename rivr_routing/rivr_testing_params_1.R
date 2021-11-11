################################################################################
### Parameters to Use ###
################################################################################

## channel slope ( vertical m / horizontal m) ##
#slopes <- 2*10^seq(-3,-2,0.25)
slopes <- 0.001
## channel length (m) ##
extents <- 5000
## pde ##
pdes <- "Dynamic" #c("Kinematic","Dynamic")
## max spatial resolution (m) ##
max_dxs <- 25
## Manning's roughness ##
mannings <- 0.1
## channel bottom width (m) ## 
widths <- 40
## channel side slope (horizontal m / vertical m) ##
sideslopes <- 1
## time step (s) ##
dts <- 5
## max time step (s) ##
max_times <- 25000
## initial flow condition (m^3/s) ##
initflows <- 25.48516
## conversion factor for Manning's equation ##
Cm <- 1.00
## gravitational acceleration (m/s^2) ##
g <- 9.81
## Set desired courant number for stability ##
# should be less than 0.06 for Dynamic wave
desired_courant_number <- 0.05
## output file template name ##
output_file_name <- '/data/rivr_outputs/outputs_testing/rivr_output_params.h5'
## aggregate files ##
aggregate_files <- FALSE
## cores to use ##
cores <- 5
## boundary condition file ##
boundary_condition_file <- '/data/usgs_data/03602500_piney_river_vernon_tn_1989_to_2021.h5'
