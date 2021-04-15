library(rivr)

### RIVR DUMMY
slope <- 0.001       # channel slope (vertical ft / horizontal ft)
extent <- 150000     # channel length (ft)
mannings <- 0.11    # Manning's roughness
width <- 100         # channel bottom width
sideslope <- 0       # channel side slope (horizontal ft / vertical ft)
Cm <- 1.486          # conversion factor for Manning's equation
g <- 32.2            # gravitational acceleration (ft/s^2)
numnodes <- 301                    # number of finite-difference nodes
dx <- extent/(numnodes - 1)        # distance between nodes (ft)
dt <- 20                          # time interval (s)

times <- seq(0, 76000, by = dt)
boundary <- ifelse(times < 9000,                  # upstream hyrograph (ft^3/s)
                   250 + (750/pi) * (1 - cos(pi * times/(4500))), 250)

monpoints <- seq(1, numnodes,1)                  # Nodes to monitor
montimes <- seq(1, length(boundary), by = 10)     # time steps to monitor
initflow <- 250                                   # initial flow condition (ft^3/s)

w <- route_wave(slope, mannings, Cm, g, width, sideslope, initflow, boundary,
                timestep = dt, spacestep = dx, numnodes = numnodes,
                monitor.nodes = monpoints, monitor.times = montimes,engine = "Kinematic")

downstream <- rep(-1, length(boundary)) 

d <- route_wave(slope, mannings, Cm, g, width, sideslope, initflow, boundary,
                downstream, timestep = dt, spacestep = dx, numnodes = numnodes,
                monitor.nodes = monpoints, monitor.times = montimes, engine = "Dynamic",
                boundary.type = "QQ")
