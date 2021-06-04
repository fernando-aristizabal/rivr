##########################################
##########################################
#### River Routing Toy Training Data #####
## Fernando Aristizabal
## In support of dynamic wave project
## Collaborating with Greg Petrochenkov

##########################################
##########################################

library(rivr)
library(doParallel)
library(foreach)
library(rhdf5)
library(stringr)

########## TO DO ##########################
# 1) Get stability
# 2) Realistic values
# 3) run them all
###########################################

## clear workspace and garbage collect
rm(list = ls())
gc()

# source input parameters
source("/rivr/rivr_testing_params_1.R")

################################################################################
################################################################################

input_params <- expand.grid( slopes, extents, mannings, widths, sideslopes, dts, max_times, initflows)

# column names
names(input_params) <- c('slope','extent_m','mannings','width_m','sideslope',
                         'dt_s','maxtime_s','initflow_cms')

# number of experiments
number_of_experiments <- dim(input_params)[1]

# add experiment number column
input_params$experiment_number <- 1:number_of_experiments

# output file parameters
ext <- tools::file_ext(output_file_name)
base_file <- tools::file_path_sans_ext(output_file_name)
id_len <- nchar(toString(number_of_experiments))
outputs_directory <- dirname(output_file_name)
if ( !dir.exists(outputs_directory) ) { dir.create(outputs_directory,recursive=TRUE, mode='00774') } 

# write input parameters to hdf5 file
if (file.exists(output_file_name)) { file.remove(output_file_name) }
h5createFile(output_file_name)
h5createGroup(output_file_name,'rivr_data')
h5createGroup(output_file_name,'rivr_data/outputs')
h5write(input_params,file=output_file_name,name='rivr_data/inputs')


################################################################################
################################################################################

# detecting number of cores on machine
#cores <- detectCores()[1]

# picks minimum of cores-1 or number of parameters
cl <- makeCluster( min( 
                       c(cores, number_of_experiments)
                      ),
                   outfile=""
                 )  

# print number of experiments
message(paste("Running",number_of_experiments,"experiments."))

# register number of cores
registerDoParallel(cl)

finalDF <- foreach (i=1:number_of_experiments) %dopar%{
#for (i in 1:number_of_experiments) {
  message(paste('Processing',i,'of',number_of_experiments))
  
  # get input parameter #
  max_time <- input_params[i,'maxtime_s']
  dt <- input_params[i,'dt_s']
  slope <- input_params[i,'slope']
  manning <- input_params[i,'mannings']
  width <- input_params[i,'width_m']
  sideslope <- input_params[i,'sideslope']
  extent <- input_params[i,'extent_m']
  initflow <- input_params[i,'initflow_cms']

  # making times and upstream boundary vectors
  times <- seq(0, max_time, by = dt)
  num_of_times <- length(times)

  # upstream hyrograph (m^3/s)
  boundary <- ifelse(times < 9000,
                     initflow + (750/pi) * (1 - cos(pi * times/(4500))), initflow)
  
  # set downstream condition (<0 is zero gradient condition)
  downstream <- rep(-1, length(boundary)) 
  
  # trying to find optimal spatial resolution, dx, to optimize courant number
  dxs <- rep(NA,length(times))
  iii <- 1
  for (b in boundary) {
    norm_dep <- rivr::normal_depth(slope, manning, b, 1, Cm, width, sideslope)
    channel_geometry <- rivr::channel_geom(norm_dep,width,sideslope)
    velocity <- b / channel_geometry[['A']]
    dxs[iii] <- velocity * dt / desired_courant_number
    iii <- iii + 1
  }
  dx <- max(dxs)
  numnodes <- ceiling(extent/dx) + 1
  dx <- extent / (numnodes - 1)
  
  # set monitoring nodes and times 
  monpoints <- 1:numnodes  # Nodes to monitor
  montimes <- 1:num_of_times  # time steps to monitor

  # solve for wave
  d <- rivr::route_wave(slope, manning, Cm, g, width, sideslope, 
                        initflow, boundary,
                        downstream, timestep = dt, spacestep = dx, numnodes = numnodes,
                        monitor.nodes = monpoints, monitor.times = montimes, engine = "Dynamic",
                        boundary.type = "QQ"
                        )
  
  # crop d data frame. for some reason it duplicates data. 
  # maybe bc of steps/nodes then monnodes/montimes
  number_of_simulation_points <- num_of_times * numnodes
  
  # crop d
  d$node[1:number_of_simulation_points]
  d$step[1:number_of_simulation_points]
  d$distance[1:number_of_simulation_points]
  d$time[1:number_of_simulation_points]
  d$flow[1:number_of_simulation_points]
  d$velocity[1:number_of_simulation_points]
  d$depth[1:number_of_simulation_points]
  d$area[1:number_of_simulation_points]
  
  # data frame building of data
  courant_numbers <- d$velocity * dt / dx
  any_na_flows <- rep(any(is.na(d$flow)),number_of_simulation_points) # if any NAs
  
  # repeats upstream and downstream boundary conditions every cycle of time indices
  all_boundaries <- rep(NA,number_of_simulation_points)
  all_downstreams <- rep(NA,number_of_simulation_points)
  for (ii in 1:number_of_simulation_points){
    s <- d$step[ii] ; n <- d$node[ii]
    idx <- s + (n-1)*num_of_times
    all_boundaries[idx] <- boundary[s]
    all_downstreams[idx] <- downstream[s]
  }
  
  # data frame
  outputs <- data.frame(
                    cbind(
                          rep(i,number_of_simulation_points),
                          any_na_flows,
                          all_boundaries,
                          all_downstreams,
                          courant_numbers[1:number_of_simulation_points],
                          d$node[1:number_of_simulation_points],
                          d$step[1:number_of_simulation_points],
                          d$distance[1:number_of_simulation_points],
                          d$time[1:number_of_simulation_points],
                          d$flow[1:number_of_simulation_points],
                          d$velocity[1:number_of_simulation_points],
                          d$depth[1:number_of_simulation_points],
                          d$area[1:number_of_simulation_points]
    )
  )
  
  # column names
  names(outputs) <- c('experiment_number',
                      'any_NAs_in_experiment','upstream_boundary_cms',
                      'downstream_m','courant_number','node','step',
                      'distance_m','time_s','flow_cms','velocity_mps',
                      'depth_m','area_sqm')
  
  # write out outputs
  id <- toString(i)
  full_id <- stringr::str_pad(id,id_len,"0",side='left')
  write.csv(outputs,paste0(base_file,"_",full_id,".csv"),row.names=F)
  
}  

# stop parallel compute cluster
stopCluster(cl)

#### aggregate CSV files to hdf5 ###
message("Aggregating CSV files to HDF5 ...")
csv_list <- Sys.glob(file.path(outputs_directory,"*.csv"))
pb <- txtProgressBar(min = 0, max= length(csv_list), style=3, char="=")

i <- 1
for (csv in csv_list) {
  outputs <- read.csv(csv)
  full_id <- tail(stringr::str_split(tools::file_path_sans_ext(basename(csv)),'_')[[1]],n=1)
  output_dataset_name <- paste0('rivr_data/outputs/',full_id)
  rhdf5::h5write(outputs,file=output_file_name,name=output_dataset_name)
  file.remove(csv)
  setTxtProgressBar(pb,i)
  i <- i + 1
}

close(pb)
message('\n')
