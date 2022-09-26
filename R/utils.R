getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

cluster_dist <- function(latent, data, start.idx) {
  cl <- data$cl
  start.clus = getmode(cl[start.idx])

  latent <- cbind(latent, cl)
  latent <- as.data.frame(latent)
  cell_stages <- factor(cl)

  x <- latent
  total_col_num <- ncol(x)
  # Euclidian center
  clus_cent <- NULL
  # Select the last column, cell type/clustering result, which has to be represent by numeric id.
  for (i in 1:length(unique(x[, total_col_num]))) {
    temp <- x[which(x[, total_col_num] == i), ]
    temp_cent <- apply(temp, 2, mean)
    clus_cent <- rbind(clus_cent, temp_cent)
  }

  clus_cent_row_names <- paste("clus", seq(1:length(unique(x[, total_col_num]))))
  clus_cent <- rbind(clus_cent[start.clus, ], clus_cent[-start.clus, ])
  clus_cent_row_names <- c(clus_cent_row_names[start.clus], clus_cent_row_names[-start.clus])
  rownames(clus_cent) <- clus_cent_row_names

  clus_dist <- as.matrix(stats::dist(clus_cent[, 1:(total_col_num - 1)], method = "euclidian"))[, 1]
  clus_cent <- as.data.frame(clus_cent)

  dist_order <- NULL
  for (i in 1:length(latent[, total_col_num])) {
    dist_order <- c(dist_order, clus_dist[which(clus_cent$cl == latent[, total_col_num][i])])
  }

  out <- list(center = clus_cent, clus_dist = clus_dist, dist_order = dist_order)
}

sort_igraph_object <- function(g, dist_vector, start_point) {
  next_nodes <- g %>%
    igraph::neighbors(start_point) %>%
    names()

  if (length(next_nodes) == 0) {
    return(list(g = g, s = start_point))
  } else {
    g_vertices <- igraph::V(g) %>% names()
    if (dist_vector[start_point] != min(dist_vector[g_vertices])) {
      min_node <- names(dist_vector)[dist_vector == min(dist_vector[g_vertices])]

      bool_adjancent <- igraph::are_adjacent(g, min_node, start_point)
      if (bool_adjancent && length(g_vertices) == 2) {
        start_point <- min_node
      } else {
        neighbour_nodes_start <- g %>%
          igraph::neighbors(start_point) %>%
          names()
        neighbour_nodes_start <- neighbour_nodes_start[!neighbour_nodes_start == min_node]
        neighbour_nodes_min <- g %>%
          igraph::neighbors(min_node) %>%
          names()
        neighbour_nodes_min <- neighbour_nodes_min[!neighbour_nodes_min == start_point]

        g <- g %>% igraph::delete_vertices(c(start_point, min_node))

        edges_start <- start_point %>%
          rep(length(neighbour_nodes_min)) %>%
          rbind(neighbour_nodes_min) %>%
          c()
        edges_min <- min_node %>%
          rep(length(neighbour_nodes_start)) %>%
          rbind(neighbour_nodes_start) %>%
          c()
        if (bool_adjancent) {
          edges_start <- edges_start %>% c(c(start_point, min_node))
        }
        g <- g %>%
          igraph::add_vertices(2, name = c(start_point, min_node)) %>%
          igraph::add_edges(edges_start) %>%
          # add_vertices(1, name = edges_min) %>%
          igraph::add_edges(edges_min)
        start_point <- min_node
      }
    }
    next_nodes <- g %>%
      igraph::neighbors(start_point) %>%
      names()
    next_g <- g %>%
      igraph::delete_vertices(start_point) %>%
      igraph::decompose.graph()

    sorted_sub_graph <- sort_igraph_object(next_g[[1]], dist_vector, intersect(next_nodes, igraph::V(next_g[[1]]) %>% names()))
    sorted_g <- sorted_sub_graph$g %>%
      igraph::add_vertices(1, name = start_point) %>%
      igraph::add_edges(c(start_point, sorted_sub_graph$s))
    if (length(next_g) > 1) {
      for (id in 2:length(next_g)) {
        temp_sub_graph <- sort_igraph_object(next_g[[id]], dist_vector, intersect(next_nodes, igraph::V(next_g[[id]]) %>% names()))
        temp_sub_graph$g <- temp_sub_graph$g %>%
          igraph::add_vertices(1, name = start_point) %>%
          igraph::add_edges(c(start_point, temp_sub_graph$s))
        sorted_g <- igraph::union(sorted_g, temp_sub_graph$g)
      }
    }
  }
  return(list(g = sorted_g, s = start_point))
}
