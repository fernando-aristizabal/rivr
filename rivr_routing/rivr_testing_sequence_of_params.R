require(rivr)

### RIVR DUMMY
num_of_params <- 2
slopes <- c(0.001,0.01)       # channel slope (vertical ft / horizontal ft)
extents <- c(150000,150000)     # channel length (ft)
mannings <- c(0.11,0.07)    # Manning's roughness
widths <- c(100,200)         # channel bottom width
sideslopes <- c(0,1)       # channel side slope (horizontal ft / vertical ft)
Cms <- c(1.486,1.486)          # conversion factor for Manning's equation
gs <- c(32.2,32.2)            # gravitational acceleration (ft/s^2)
#numnodes <- c(301,151)                    # number of finite-difference nodes
#dxs <- extents/(numnodes - 1)        # distance between nodes (ft)
dts <- c(20,20)                          # time interval (s)
max_times <- c(76000,76000)
initflows <- c(250,250)                                   # initial flow condition (ft^3/s)

desired_courant_number <- 0.04 # should be less than 0.06
#normal_depth <- normal_depth(slopes[2],mannings[2],initflows[2],1,Cms[2],widths[2],sideslopes[2])
#channel_geometry <- channel_geom(normal_depth,widths[2],sideslopes[2])
#velocity <- initflows[2] / channel_geometry[['A']]

#courant_number <- velocity * dts[2] / dxs[2]
#dx <- velocity * dts[2] / max_courant_number

for (i in 1:num_of_params) {
  print(paste('Processing',i,'of',num_of_params))
  times <- seq(0, max_times[i], by = dts[i])
  boundary <- ifelse(times < 9000,                  # upstream hyrograph (ft^3/s)
                     250 + (750/pi) * (1 - cos(pi * times/(4500))), 250)
  
  dxs <- rep(NA,length(times))
  iii <- 1
  for (b in boundary) {
    normal_depth <- normal_depth(slopes[i],mannings[i],b,1,Cms[2],widths[i],sideslopes[i])
    channel_geometry <- channel_geom(normal_depth,widths[i],sideslopes[i])
    velocity <- b / channel_geometry[['A']]
    dxs[iii] <- velocity * dts[i] / desired_courant_number
    iii <- iii + 1
  }

  dx <- max(dxs)
  numnodes <- ceiling(extents[i]/dx) + 1
  dx <- extents[i]/(numnodes - 1)
  
  monpoints <- seq(1, numnodes,1)                  # Nodes to monitor
  montimes <- seq(1, length(boundary), by = dts[i])     # time steps to monitor

  downstream <- rep(-1, length(boundary)) # downstream condition (<0 is zero gradient condition)

  print('   Routing wave')
  d <- route_wave(slopes[i], mannings[i], Cms[i], gs[i], widths[i], sideslopes[i], 
                  initflows[i], boundary,
                  downstream, timestep = dts[i], spacestep = dx, numnodes = numnodes,
                  monitor.nodes = monpoints, monitor.times = montimes, engine = "Dynamic",
                  boundary.type = "QQ")
  
  # data frame
  print('   Building dataframe')
  num_of_sims <- length(d$flow)
  num_of_times <- length(times)
  all_boundaries <- rep(NA,num_of_sims)
  all_downstreams <- rep(NA,num_of_sims)
  
  courant_numbers <- d$velocity * dts[i] / dx

  for (ii in 1:num_of_sims){
    s <- d$step[ii] ; n <- d$node[ii]
    idx <- s + (n-1)*num_of_times
    all_boundaries[idx] <- boundary[s]
    all_downstreams[idx] <- downstream[s]
  }

  df <- data.frame(
                    cbind(rep(slopes[i],num_of_sims),
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
  
  names(df) <- c('slope','mannings','width_ft','initflow_cfs',
                 'dx_ft','dt_s','upstream_boundary_cfs','downstream_ft',
                 'courant_number','node','step','distance_ft',
                 'time_s','flow_cfs','velocity_fps','depth_ft',
                 'area_sqft')
  
  if (i==1){
    dfs <- df
  } else {
    dfs <- rbind(dfs,df)
  }
  
}  


