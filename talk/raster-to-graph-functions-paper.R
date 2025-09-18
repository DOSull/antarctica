library(sf)
library(terra)
library(igraph)
library(tidyr)
library(dplyr)

# The tobler hiking function for speed appropriately scaled to m/hr
# if modified, use Marquez-Perez modified form which is a bit less
# superman-ish. These formulae from movecost help. See also 
# 
# Tobler WR. 1993. Three Presentations on Geographical Analysis and Modeling: 
#   Non-Isotropic Geographic Modeling; Speculations on the Geometry of Geography;
#   and Global Spatial Analysis. National Center for Geographic Information and 
#   Analysis, Santa Barbara, CA. https://escholarship.org/uc/item/05r820mz
#
tobler_speed_m_hr <- function(rise_over_run, modified = TRUE) {
  if (modified) {
    # Tobler off-path
    3600 * exp(-3.5 * abs(rise_over_run + 0.05))
  } else {
    # Tobler - superman!
    6 * exp(-3.5 * abs(rise_over_run + 0.05)) * 1000
  }
}

est_hiking_function <- function(slope, terrain1, terrain2) {
  if (terrain1 == terrain2) {
    if (terrain1 == "moraine") {
      4169.969 * exp(-(slope + 0.05259) ^ 2 / 0.235957)
    } else {
      3763.969 * exp(-(slope + 0.1188) ^ 2 / 0.364931)
    }
  } else {
    (4169.969 * exp(-(slope + 0.05259) ^ 2 / 0.235957) + 3763.969 * exp(-(slope + 0.1188) ^ 2 / 0.364931)) / 2
  }
}

est_hiking_function_vec <- Vectorize(est_hiking_function)

est_hiking_function_z <- function(slope, terrain1, terrain2) {
  if (terrain1 == terrain2) {
    if (terrain1 == "moraine") {
      4235.085 * exp(-(slope + 0.061365) ^ 2 / 0.1266212)
    } else {
      3870.854 * exp(-(slope + 0.118810) ^ 2 / 0.2975753)
    }
  } else {
    (4235.085 * exp(-(slope + 0.061365) ^ 2 / 0.1266212) + 3870.854 * exp(-(slope + 0.118810) ^ 2 / 0.2975753)) / 2
  }
}

est_hiking_function_z_vec <- Vectorize(est_hiking_function_z)

# Get a hillshade raster from terrain
get_hillshade <- function(r) {
  aspect <- terrain(r, "aspect", unit = "radians")
  slope <- terrain(r, "slope", unit = "radians")
  shade(slope, aspect, 40, 135)
}

# convert a terra raster to a sf points dataset
raster_to_points <- function(raster) {
  raster |>
    as.points() |>
    st_as_sf() |>
    st_sf()
}

# # get x, y, z coordinates from a sf points dataset
# # there is assumed to be one additional variate, the elevation in column 1
# extract_xyz_from_points <- function(pts) {
#   xy <- st_coordinates(pts)
#   list(x = xy[, "X"],
#        y = xy[, "Y"],
#        z = st_drop_geometry(pts)[, 1])
# }

# get an igraph graph from points using distance as the adjacency criterion
graph_from_points <- function(pts, distance) {
  pts |>
    # st_intersects on a buffered copy is MUCH quicker than st_is_within_distance
    st_intersects(pts |> st_buffer(distance)) |>
    # going via sparse matrix required to avoid memory running out on full matrices
    as("sparseMatrix") |>
    # no self-loops
    graph_from_adjacency_matrix(diag = FALSE)
}

# make a list with $i and $j members for from-to indices of edges of graph G 
get_graph_adjacencies <- function(G) {
  edges <- ends(G, E(G))
  list(i = edges[, 1], j = edges[, 2])
}

# assign a range of useful attributes to the nodes and edges of graph G
# takes a vector of movement costs (impedances) in same order as list xyz
# with x, y, z members for the spatial coordinates
assign_movement_variables_to_graph_2 <- function(G, xyz, use_gps_elevations = TRUE,
                                                 terrain = rep("terrain", length(G)),
                                                 impedances = rep(1, length(G))) {
  V(G)$x <- xyz$x
  V(G)$y <- xyz$y
  V(G)$z <- xyz$z
  # get graph edge indice to index into x, y, z for the calculations
  ij <- get_graph_adjacencies(G)
  i <- ij$i  # the 'from' index of each edge
  j <- ij$j  # the 'to' index of each edge
  E(G)$length_xy <- sqrt((V(G)$x[i] - V(G)$x[j]) ^ 2 + 
                         (V(G)$y[i] - V(G)$y[j]) ^ 2)
  E(G)$z_diff <- V(G)$z[i] - V(G)$z[j]
  E(G)$gradient <- E(G)$z_diff / E(G)$length_xy
  if (use_gps_elevations) {
    E(G)$speed <- est_hiking_function_vec(E(G)$gradient, terrain[i], terrain[j])
  } else {
    E(G)$speed <- est_hiking_function_z_vec(E(G)$gradient, terrain[i], terrain[j])
  }
  E(G)$cost <- E(G)$length_xy / E(G)$speed
  G
}

get_rgb_from_imagery <- function(image_raster, points) {
  rgb <- terra::extract(image_raster, as(points, "SpatVector")) 
  names(rgb)[2:4] <- c("R", "G", "B")
  rgb
}

get_xy1_to_xy2_transform_df <- function(G, origin_index, rgb) {
  x0 <- V(G)$x[origin_index]
  y0 <- V(G)$y[origin_index]
  data.frame(
    x1 = V(G)$x, y1 = V(G)$y, z = V(G)$z,
    red = rgb$R, green = rgb$G, blue = rgb$B) |>
    mutate(
      ID = row_number(), dx = x1 - x0, dy = y1 - y0,
      x1 = round(x1, 1), y1 = round(y1, 1), z = round(z, 1),
      distance = round(sqrt(dx ^2 + dy ^ 2), 1),
      bearing = atan2(dy, dx),
      time_distance = round(V(G)$time_hrs * 900, 1),
      x2 = round(x1 + time_distance * cos(bearing), 1), 
      y2 = round(y1 + time_distance * sin(bearing), 1)) |>
    # drop dx, dy, and bearing
    select(-dx, -dy, -bearing) |>
    # reorder by time_distance and assign an index ordering
    arrange(time_distance) |>
    mutate(z_order = row_number()
    )
}

make_linestring <- function(x1, y1, x2, y2) {
  st_linestring(matrix(c(x1, y1, x2, y2), ncol = 2, byrow = TRUE))
}

get_graph_as_line_layer <- function(G, crs = 3031, append_attributes = FALSE) {
  ij <- get_graph_adjacencies(G)
  i <- ij$i
  j <- ij$j
  edges_as_sf <- mapply(make_linestring, 
                        V(G)$x[i], V(G)$y[i], V(G)$x[j], V(G)$y[j],
                        SIMPLIFY = FALSE) |>
    st_sfc() |>
    st_as_sf(crs = crs) |>
    rename(geometry = x)
  if (length(edge_attr_names(G)) > 0 & append_attributes) {
    for (attr in edge_attr_names(G)) {
      edges_as_sf[attr] <- edge_attr(G, attr)
    }
  }
  edges_as_sf                      
}

get_graph_as_point_layer <- function(G, crs = 3031) {
  ij <- get_graph_adjacencies(G)
  i <- ij$i
  j <- ij$j
  vertices_as_sf <- matrix(c(V(G)$x, V(G)$y, V(G)$z), ncol = 3) |>
    st_multipoint() |>
    st_sfc() |>
    st_cast("POINT") |>
    st_as_sf(crs = crs) |>
    rename(geometry = x)
  if (length(vertex_attr_names(G)) > 0) {
    for (attr in vertex_attr_names(G)) {
      vertices_as_sf[attr] <- vertex_attr(G, attr)
    }
  }
  vertices_as_sf                      
}

get_shortest_path_tree <- function(G, vid) {
  j <- V(G) |> 
    as.vector()
  i <- shortest_paths(
    G, from = vid, to = j, mode = "out", weights = edge_attr(G, "cost"),
    predecessors = TRUE)$predecessors |>
    as.vector()
  i <- i[-vid]
  j <- j[-vid]
  graph_from_edgelist(matrix(c(i, j), ncol = 2)) |>
    set_vertex_attr("x", value = vertex_attr(G, "x")) |>
    set_vertex_attr("y", value = vertex_attr(G, "y"))
}

# see https://rmazing.wordpress.com/2013/08/14/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/
predictNLS <- function(object, newdata, level = 0.95, nsim = 10000, ...) {
  require(MASS, quietly = TRUE)
  ## get right-hand side of formula
  RHS <- as.list(object$call$formula)[[3]]
  EXPR <- as.expression(RHS)
  ## all variables in model
  VARS <- all.vars(EXPR)
  ## coefficients
  COEF <- coef(object)
  ## extract predictor variable    
  predNAME <- setdiff(VARS, names(COEF))  
  ## take fitted values, if 'newdata' is missing
  if (missing(newdata)) {
    newdata <- eval(object$data)[predNAME]
    colnames(newdata) <- predNAME
  }
  ## check that 'newdata' has same name as predVAR
  if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")
  ## get parameter coefficients
  COEF <- coef(object)
  ## get variance-covariance matrix
  VCOV <- vcov(object)
  ## augment variance-covariance matrix for 'mvrnorm' 
  ## by adding a column/row for 'error in x'
  NCOL <- ncol(VCOV)
  ADD1 <- c(rep(0, NCOL))
  ADD1 <- matrix(ADD1, ncol = 1)
  colnames(ADD1) <- predNAME
  VCOV <- cbind(VCOV, ADD1)
  ADD2 <- c(rep(0, NCOL + 1))
  ADD2 <- matrix(ADD2, nrow = 1)
  rownames(ADD2) <- predNAME
  VCOV <- rbind(VCOV, ADD2) 
  ## iterate over all entries in 'newdata' as in usual 'predict.' functions
  NR <- nrow(newdata)
  respVEC <- numeric(NR)
  seVEC <- numeric(NR)
  varPLACE <- ncol(VCOV)   
  ## define counter function
  counter <- function (i) {
    if (i%%10 == 0) 
      cat(i)
    else cat(".")
    if (i%%50 == 0) 
      cat("\n")
    flush.console()
  }
  outMAT <- NULL 
  for (i in 1:NR) {
    counter(i)
    ## get predictor values and optional errors
    predVAL <- newdata[i, 1]
    if (ncol(newdata) == 2) predERROR <- newdata[i, 2] else predERROR <- 0
    names(predVAL) <- predNAME  
    names(predERROR) <- predNAME  
    ## create mean vector for 'mvrnorm'
    MU <- c(COEF, predVAL)
    ## create variance-covariance matrix for 'mvrnorm'
    ## by putting error^2 in lower-right position of VCOV
    newVCOV <- VCOV
    newVCOV[varPLACE, varPLACE] <- predERROR^2
    ## create MC simulation matrix
    simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
    ## evaluate expression on rows of simMAT
    EVAL <- try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
    if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
    ## collect statistics
    PRED <- data.frame(predVAL)
    colnames(PRED) <- predNAME   
    FITTED <- predict(object, newdata = data.frame(PRED))
    MEAN.sim <- mean(EVAL, na.rm = TRUE)
    SD.sim <- sd(EVAL, na.rm = TRUE)
    MEDIAN.sim <- median(EVAL, na.rm = TRUE)
    MAD.sim <- mad(EVAL, na.rm = TRUE)
    QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
    RES <- c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
    outMAT <- rbind(outMAT, RES)
  }
  colnames(outMAT) <- c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
  rownames(outMAT) <- NULL
  cat("\n")
  return(outMAT)  
}
