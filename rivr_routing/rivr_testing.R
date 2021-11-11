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

###########################################
###########################################

## clear workspace and garbage collect
rm(list = ls())
gc()

# source input parameters
source("/rivr/rivr_testing_params_3.R")

# source functions
source("rivr/rivr_functions.R")

# Set UMASK
Sys.umask(mode="002")

################################################################################
################################################################################

experiment_params <- expand.grid( slopes, extents, max_dxs, mannings, widths, sideslopes, dts, max_times, initflows, pdes, stringsAsFactors=FALSE)

# column names
names(experiment_params) <- c('slope','extent_m','dx_m','mannings','width_m','sideslope',
                              'dt_s','maxtime_s','initflow_cms','pde')

# load boundary condition if using empirical datasource
if (!is.na(boundary_condition_file)) {
  boundary <- h5read(boundary_condition_file,'usgs_data_processed/table',bit64conversion='double')
  datetimes <- boundary[['datetime']]
  boundary <- boundary[['streamflow_cms']]
  datetimes <- strptime(datetimes,format='%Y-%m-%d %H:%M:%S',tz='CST')
  
  #maxtime_s <- as.double(difftime(datetimes[length(datetimes)],datetimes[1],units='secs'))
  #dt_s <- as.double(difftime(datetimes[2],datetimes[1],units='secs'))
  maxtime_s <- length(datetimes) - 1
  dt_s <- 1

  experiment_params$maxtime_s <- maxtime_s
  experiment_params$dt_s <- dt_s

  experiment_params <- experiment_params[!duplicated(experiment_params),]
  
}

# number of experiments
number_of_experiments <- dim(experiment_params)[1]

# add experiment number column
experiment_params$experiment_number <- 1:number_of_experiments

# output file parameters
ext <- tools::file_ext(output_file_name)
base_file <- tools::file_path_sans_ext(output_file_name)
id_len <- nchar(toString(number_of_experiments))
outputs_directory <- dirname(output_file_name)
if ( !dir.exists(outputs_directory) ) { dir.create(outputs_directory,recursive=TRUE, mode='00774') } 

# create hdf5 file and groups
if (file.exists(output_file_name)) { file.remove(output_file_name) }
h5createFile(output_file_name)
h5createGroup(output_file_name,'rivr_data')
h5createGroup(output_file_name,'rivr_data/outputs')


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

additional_experiment_params <- foreach (i=1:number_of_experiments,.combine='rbind') %dopar% {
#for (i in 1:number_of_experiments) {
    
    # message update
    message(paste('Processing',i,'of',number_of_experiments))

    # get input parameter #
    max_time <- experiment_params[i,'maxtime_s']
    dt <- experiment_params[i,'dt_s']
    slope <- experiment_params[i,'slope']
    manning <- experiment_params[i,'mannings']
    pde <- experiment_params[i,'pde']
    width <- experiment_params[i,'width_m']
    sideslope <- experiment_params[i,'sideslope']
    extent <- experiment_params[i,'extent_m']
    max_dx <- experiment_params[i,'dx_m']
    initflow <- experiment_params[i,'initflow_cms']

    # making times and upstream boundary vectors
    times <- seq(0, max_time, by = dt)
    num_of_times <- length(times)

    # set downstream condition (<0 is zero gradient condition)
    downstream <- rep(-1, length(times)) 

    # upstream hyrograph (m^3/s)
    #boundary <- ifelse(((times >= 5000) & (times < 10000)) | ((times>=15000) & (times<20000)),
    #                   initflow + (750/pi) * (1 - cos(pi * times/(4500))), initflow)

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

    dx <- max_dx

    numnodes <- ceiling(extent/dx) + 1
    dx <- extent / (numnodes - 1)
    #print(paste(dx,numnodes,dt,num_of_times))

    # set monitoring nodes and times 
    monpoints <- 1:numnodes  # Nodes to monitor
    montimes <- 1:num_of_times  # time steps to monitor

    # solve for wave and time
    start_time <- Sys.time()
    d <- rivr::route_wave(slope, manning, Cm, g, width, sideslope, 
                        initflow, boundary,
                        downstream, timestep = dt, spacestep = dx, numnodes = numnodes,
                        monitor.nodes = monpoints, monitor.times = montimes, engine = pde,
                        boundary.type = "QQ"
                        )
    compute_time_s <- as.double(Sys.time() - start_time,units='secs')

    # crop d data frame. for some reason it duplicates data. 
    # maybe bc of steps/nodes then monnodes/montimes
    number_of_simulation_points <- num_of_times * numnodes

    # crop d
    d <- d[1:number_of_simulation_points,]
    #d$node[1:number_of_simulation_points]
    #d$step[1:number_of_simulation_points]
    #d$distance[1:number_of_simulation_points]
    #d$time[1:number_of_simulation_points]
    #d$flow[1:number_of_simulation_points]
    #d$velocity[1:number_of_simulation_points]
    #d$depth[1:number_of_simulation_points]
    #d$area[1:number_of_simulation_points]

    # determine values: courant numbers & friction slope
    courant_numbers <- d$velocity * dt / dx
    #hydraulic_radius <- sapply(lapply(d$depth,rivr::channel_geom,B=width,SS=sideslope),function(x) as.numeric(x['R']))
    hydraulic_radius <- sapply(d$depth,function(D,B,SS) as.numeric(rivr::channel_geom(D,B,SS)['R']),B=width,SS=sideslope)
    friction_slope <- friction_slope_mannings(mannings,d$velocity,Cm,hydraulic_radius)
    #print(max(courant_numbers[!is.na(courant_numbers)]))

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
                         experiment_number = rep(i,number_of_simulation_points),
                         upstream_boundary_cms = all_boundaries,
                         downstream_m = all_downstreams,
                         courant_number = courant_numbers[1:number_of_simulation_points],
                         node = d$node[1:number_of_simulation_points],
                         step = d$step[1:number_of_simulation_points],
                         distance_m = d$distance[1:number_of_simulation_points],
                         time_s = d$time[1:number_of_simulation_points],
                         flow_cms = d$flow[1:number_of_simulation_points],
                         velocity_mps = d$velocity[1:number_of_simulation_points],
                         depth_m = d$depth[1:number_of_simulation_points],
                         area_sqm = d$area[1:number_of_simulation_points],
                         friction_slope = friction_slope,
                         hydraulic_radius = hydraulic_radius
                        )
    #print(sum(is.na(outputs$flow_cms))/length(outputs$flow_cms))

    # experiment ids
    id <- toString(i)
    full_id <- stringr::str_pad(id,id_len,"0",side='left')

    # create hdf5 file and groups
    experiment_output_file_name <- paste0(base_file,"_",full_id,".h5")
    if (file.exists(experiment_output_file_name)) { file.remove(experiment_output_file_name) }
    rhdf5::h5write(outputs,file=experiment_output_file_name,name='rivr_data')

    # plotting boundary conditions
    #jpeg(file=paste0("/data/boundary_",full_id,".jpg"))
    #plot(times,log10(boundary),xlab='Times (s)',ylab='Discharge (CMS)',main=paste('Boundary Condition for Experiment',full_id),type='l',log='y')
    #dev.off()

    # get stability of experiment by checking for nas
    stability <- as.integer(!any(is.na(d$flow)))
    message(paste("Exp#:",i ,"| Stable:", stability, "| Compute Time:", compute_time_s,"| dx:",dx,
                "| mannings:",mannings, "| slope:", slope, "| sideslope:",sideslope,"| width:",width))
    message('')

    # create dataframe with additional experiment parameters to rbind 
    #additional_experiment_params <- data.frame( 
    data.frame( 
              stability=stability,
              dx_m=dx,
              num_of_nodes=numnodes,
              num_of_times = num_of_times,
              compute_time_s = compute_time_s
            )

}  

# stop parallel compute cluster
stopCluster(cl)

# write experiment params
experiment_params <- subset(experiment_params,select=-dx_m)
experiment_params <- cbind(experiment_params,additional_experiment_params)
h5write(experiment_params,file=output_file_name,name='rivr_data/experiments')

# Summarize stability
message(paste("Stable experiments:",sum(experiment_params$stability),"out of",number_of_experiments))


#### aggregate CSV files to hdf5 ###
if (aggregate_files){
    message("Aggregating HDF5 files ...")
    experiment_output_file_names_to_glob <- paste0(base_file,"_*.h5")
    file_list <- Sys.glob(experiment_output_file_names_to_glob)
    pb <- txtProgressBar(min = 0, max= length(file_list), style=3, char="=")

    i <- 1
    for (f in file_list) {
      outputs <- h5read(f,'rivr_data')
      full_id <- tail(stringr::str_split(tools::file_path_sans_ext(basename(f)),'_')[[1]],n=1)
      output_dataset_name <- paste0('rivr_data/outputs/',full_id)
      h5write(outputs,file=output_file_name,name=output_dataset_name)
      file.remove(f)
      setTxtProgressBar(pb,i)
      i <- i + 1
    }
    
    close(pb)
}


message('')
