library(rivr)
library(doParallel)
library(foreach)


########## TO DO ##########################
# 3) Play with QQ/solver details
# 4) save engine/QQ/solver too
# 5) build larger sequences
# 6) run them all
###########################################

## clear workspace and garbage collect
rm(list = ls())
gc()

#########################################################
### Parameters to Use ###

# number of experiments
number_of_experiments <- 6

## channel slope (%, vertical m / horizontal m) ##
slopes <- c(0.1,0.2,0.4,0.6,0.8,1.0) / 100
## channel length (m) ##
extents <- rep(1500,number_of_experiments)
## Manning's roughness ##
mannings <- c(0.10,0.11,0.12,0.13,0.14,0.15)
## channel bottom width (m) ## 
widths <- c(35,35,30,30,25,25)
## channel side slope (horizontal m / vertical m) ##
sideslopes <- rep(0,number_of_experiments)
## conversion factor for Manning's equation ##
Cms <- c(1.00,number_of_experiments)
## gravitational acceleration (m/s^2) ##
gs <- rep(9.81,number_of_experiments)
## time step (s) ##
dts <- rep(20,number_of_experiments)
## max time step (s) ##
max_times <- rep(40000,number_of_experiments)
## initial flow condition (m^3/s) ##
initflows <- c(10,10,8,8,7,7)

## Set desired courant number for stability ##
# should be less than 0.06 for Dynamic wave
desired_courant_number <- 0.04 

################################################################################

# detecting number of cores on machine
cores <- detectCores()[1]

# picks minimum of cores-1 or number of parameters
cl <- makeCluster( min( 
                       c(cores[1]-1, number_of_experiments)
                      )
                 )  

# register number of cores
registerDoParallel(cl)

finalDF <- foreach (i=1:number_of_experiments, .combine=rbind) %dopar%{
  print(paste('Processing',i,'of',number_of_experiments))
  
  # making times and upstream boundary vectors
  times <- seq(0, max_times[i], by = dts[i])

  # upstream hyrograph (m^3/s)
  boundary <- ifelse(times < 9000,
                     7 + (750/pi) * (1 - cos(pi * times/(4500))), 7)
  
  # set downstream condition (<0 is zero gradient condition)
  downstream <- rep(-1, length(boundary)) 
  
  # trying to find optimal spatial resolution, dx, to optimize courant number
  dxs <- rep(NA,length(times))
  iii <- 1
  for (b in boundary) {
    norm_dep <- rivr::normal_depth(slopes[i],mannings[i],b,1,Cms[2],widths[i],sideslopes[i])
    channel_geometry <- rivr::channel_geom(norm_dep,widths[i],sideslopes[i])
    velocity <- b / channel_geometry[['A']]
    dxs[iii] <- velocity * dts[i] / desired_courant_number
    iii <- iii + 1
  }
  dx <- max(dxs)
  numnodes <- ceiling(extents[i]/dx) + 1
  dx <- extents[i]/(numnodes - 1)
 
  # set monitoring nodes and times 
  monpoints <- seq(1, numnodes,1)                    # Nodes to monitor
  montimes <- seq(1, length(boundary), by = dts[i])  # time steps to monitor

  # solve for wave
  #print('   Routing wave')
  d <- rivr::route_wave(slopes[i], mannings[i], Cms[i], gs[i], widths[i], sideslopes[i], 
                        initflows[i], boundary,
                        downstream, timestep = dts[i], spacestep = dx, numnodes = numnodes,
                        monitor.nodes = monpoints, monitor.times = montimes, engine = "Dynamic",
                        boundary.type = "QQ")
  
  # data frame building of data
  #print('   Building dataframe')
  num_of_sims <- length(d$flow)
  num_of_times <- length(times)
  courant_numbers <- d$velocity * dts[i] / dx

  # repeats upstream and downstream boundary conditions every cycle of time indices
  all_boundaries <- rep(NA,num_of_sims)
  all_downstreams <- rep(NA,num_of_sims)
  for (ii in 1:num_of_sims){
    s <- d$step[ii] ; n <- d$node[ii]
    idx <- s + (n-1)*num_of_times
    all_boundaries[idx] <- boundary[s]
    all_downstreams[idx] <- downstream[s]
  }

  # data frame
  dafr <- data.frame(
                    cbind(
                          rep(i,num_of_sims),
                          rep(slopes[i],num_of_sims),
                          rep(mannings[i],num_of_sims),
                          rep(widths[i],num_of_sims),
                          rep(initflows[i],num_of_sims),
                          rep(dxs[i],num_of_sims),
                          rep(dts[i],num_of_sims),
                          all_boundaries,
                          all_downstreams,
                          courant_numbers,
                          d$node,
                          d$step,
                          d$distance,
                          d$time,
                          d$flow,
                          d$velocity,
                          d$depth,
                          d$area
    )
  )
  
  # column names
  names(dafr) <- c('experiment_number','slope','mannings','width_m','initflow_cms',
                 'dx_m','dt_s','upstream_boundary_cms','downstream_m',
                 'courant_number','node','step','distance_m',
                 'time_s','flow_cms','velocity_mps','depth_m',
                 'area_sqm')
  dafr
  
}  

# stop parallel compute cluster
stopCluster(cl)

# check for NA flow values, presumably from unstable results
na_indices <- is.na(finalDF$flow_cms)
count_of_nas <- sum(na_indices)
total_entries <- length(finalDF$flow_cms)
message(paste(count_of_nas, 'NA flow values out',total_entries,'total flow values.'))

# write out combined df
write.csv(finalDF,'/data/rivr_params_20210428.csv',row.names=F)
