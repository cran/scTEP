#' @title preprocessing
#' @description Conduct preprocessing, including remove all zero columns and scale gene expression smaller than 100 by log transformation with 2 as base.
#' @param expr The gene expression matrix.
#' @return List with the following keys:
#' \itemize{
#' \item expr - Gene expression matrix, with rows represent samples and columns represent genes.
#' }
#' @examples
#' #Load the package
#' library(scTEP)
#' #Load example data
#' data("goolam")
#' #Get data matrix
#' expr <- as.matrix(t(SummarizedExperiment::assay(goolam)))
#'
#' data = preprocessing(expr)
#' @export
preprocessing <- function(expr) {
  expr = as.matrix(expr)
  expr <- expr[, colSums(expr) > 0]
  colnames(expr) <- toupper(colnames(expr))
  if (max(expr) > 100) expr <- log2(expr + 1)

  data <- list(
    expr = expr
  )
  return(data)
}


#' @title scTEP.fa
#' @description The 'scTEP.fa' function first selects the corresponding pathway gene sets of the data set from KEGG, then intersect the genes in the expression matrix with each pathway to have an intersect gene expression matrix for all pathways.
#' @param data A list consists of gene expression matrix.
#' @param genesets A list consists of Homo sapiens and Mus musculus gene sets.
#' @param data_org The organism of the data set, mmu or hsa.
#' @param ncores Number of processor cores to use. This values is set to \code{seed = 10L} by default
#' @param seed A parameter to set a seed for reproducibility.
#' @return List with the following keys:
#' \itemize{
#' \item faData - A large matrix consists of concatenated 2 dimensional factor analysis results of pathways.
#' }
#' @examples
#' # Load the package and the example data (goolam datas set)
#' library(scTEP)
#' #Load pathway genesets
#' data('genesets')
#' #Load example data (SCE dataset)
#' data("goolam")
#' #Get data matrix and label
#' expr <- as.matrix(t(SummarizedExperiment::assay(goolam)))[1:10, 1:100]
#'
#' #Get data matrix and label
#' data = preprocessing(expr)
#'
#' #Generate factor analysis results for all the intersections between data matrix and genesets
#' data_fa = scTEP.fa(data, genesets, ncores = 2, data_org = 'mmu', seed = 1)
#' @export
scTEP.fa <- function(data, genesets, data_org = "hsa", ncores = 10L, seed = NULL) {
  set.seed(seed)
  cl <- parallel::makeCluster(min(3, ncores), outfile = "/dev/null")
  doParallel::registerDoParallel(cl, cores = ncores)
  parallel::clusterEvalQ(cl, {
    library(psych)
  })

  datByGeneset <- lapply(genesets[[data_org]], function(genes) {
    commonGenes <- intersect(genes, colnames(data$expr))
    if (length(commonGenes) < 10) {
      return(NULL)
    }
    as.matrix(data$expr[, commonGenes])
  })

  datByGeneset <- datByGeneset[!sapply(datByGeneset, is.null)]

  dat = NULL
  . = NULL
  faData <- foreach::foreach(dat = datByGeneset) %dopar%
    {
      r <- suppressMessages(suppressWarnings(psych::fa(dat, nfactors = 2, warnings = F)$scores))
      r[r > 5] <- 5
      r[r < -5] <- -5
      r
    } %>%
    do.call(what = cbind) %>%
    `-`(min(.))
  parallel::stopCluster(cl)
  return(faData)
}



#' @importFrom dplyr %>%
#' @importFrom foreach %dopar%
#' @importFrom tibble tibble
#' @title trajectoryinference
#' @description This is the main function that performs sc-RNA seq data trajectory inference.
#' The 'scTEP' first load the latent representation and clustering result of scDHA.
#' Second, the 'trajectoryinference' function iterates through all the clustering results and calculates the distance between clusters as pseudotimes.
#' Third, it calculates the average pseudotime for every cluster in the clustering result obtained from the first step.
#' Fourth, it generates an MST and fine-tunes it by the pseudotime. Therefore, we have the trajectory for the data set.
#' @param data A list consists of gene expression matrix.
#' @param start.idx The indexes of the start cells, given by user.
#' @param scDHA_res The 'scDHA' results, consists of latent and clustering result.
#' @param allCluster A list consists of clustering results using 'scDHA' with \code{k = 5:10}.
#' @param ncores Number of processor cores to use. This values is set to \code{seed = 10L} by default.
#' @param seed A parameter to set a seed for reproducibility.
#' @return List with the following keys:
#' \itemize{
#' \item pseudotime - The pseudotime of cells in the data set.
#' \item cluster  - The clustering results of data set.
#' \item data_clus_cent - The center of all the clusters.
#' \item milestone_network - The milestone network of the inferred trajectory.
#' \item g - An 'igraph' object of the inferred trajectory
#' }
#'
#' @references
#'
#' 1. Duc Tran, Hung Nguyen, Bang Tran, Carlo La Vecchia, Hung N. Luu, Tin Nguyen (2021). Fast and precise single-cell data analysis using a hierarchical autoencoder. Nature Communications, 12, 1029. doi: 10.1038/s41467-021-21312-2
#
#' @examples
#' \donttest{
#' # Load the package and the example data (goolam data set)
#' library(scTEP)
#' #Load pathway genesets
#' data('genesets')
#' #Load example data (SCE dataset)
#' data("goolam")
#' #Get data matrix and label
#' expr <- as.matrix(t(SummarizedExperiment::assay(goolam)))
#' label <- as.character(goolam$label)
#' stages = goolam@metadata$cell.stages
#'
#' #Get data matrix and label
#' data = preprocessing(expr)
#'
#' #Generate factor analysis results for all the intersections between data matrix and genesets
#' data_fa = scTEP.fa(data, genesets, ncores = 2, data_org = 'mmu', seed = 1)
#'
#' #Get clustering results using 'scDHA' with k from 6 to 10.
#' allCluster = scTEP::clustering(data, ncores = 2)
#'
#' scDHA_res <- scDHA(data_fa, do.clus = T, gen_fil = T, ncores = 2, seed = 1)
#' #Conduct trajectory inference to the data matrix
#' idx = which(label == stages[1])
#' out = trajectoryinference(data, idx, scDHA_res, allCluster, ncores = 2, seed = 1)
#'
#' }
#' @export
trajectoryinference <- function(data, start.idx, scDHA_res, allCluster, ncores = 10L, seed = NULL) {
  set.seed(seed)

  pseudotime <- rep(0, nrow(data$expr))
  latent <- scDHA_res$latent
  df <- as.data.frame(latent)

  for (clus in allCluster) {
    data$cl <- clus
    out <- cluster_dist(latent, data, start.idx)
    tmp <- out$dist_order %>% as.numeric()
    tmp <- tmp / max(tmp)
    pseudotime <- pseudotime + tmp
  }

  # process for graph
  cl <- scDHA_res$cluster
  # Assign clusters with only 1 cells to closest group
  clus_id_remove <- which(table(cl) == 1)
  for (i in clus_id_remove) {
    temp_index <- which(cl == i)
    close_cell_index <- which.min(abs(pseudotime[-temp_index] - pseudotime[temp_index]))
    if (close_cell_index >= temp_index) {
      close_cell_index <- close_cell_index + 1
    }
    cl[temp_index] <- cl[close_cell_index]
  }

  for (i in 1:length(clus_id_remove)) {
    temp_clus_id <- clus_id_remove[i]
    cl[which(cl > temp_clus_id)] <- cl[which(cl > temp_clus_id)] - 1
    clus_id_remove <- clus_id_remove - 1
  }
  # Select the start cluster in cl
  start.clus = getmode(cl[start.idx])

  data_clus_cent <- NULL
  for (cl_id in unique(cl)) {
    cl_index <- which(cl == cl_id)
    clus_pseudotime <- mean(pseudotime[cl_index])
    tmp_clus_cent <- scDHA_res$latent[cl_index, ] %>%
      colMeans() %>%
      c(clus_pseudotime)
    data_clus_cent <- rbind(data_clus_cent, tmp_clus_cent)
  }

  data_clus_cent <- data_clus_cent %>% as.matrix()
  rownames(data_clus_cent) <- unique(cl) %>% as.character()


  distM <- data_clus_cent[, 1:15] %>%
    stats::dist() %>%
    as.matrix() # distance matrix for cluster centers

  # Build a MST tree which is consistent with pseudotime
  sorted_graph <- igraph::graph_from_adjacency_matrix(distM, weighted = TRUE, mode = "undirected") %>%
    igraph::mst() %>%
    sort_igraph_object(data_clus_cent[, 16], as.character(start.clus))

  # assign weight according to distM
  igraph::E(sorted_graph$g)$weight <- sorted_graph$g %>%
    igraph::get.edgelist() %>%
    apply(1, function(x) distM[as.numeric(x)[1], as.numeric(x)[2]])

  # Calculate the length of edge according to the euclidean distance.
  edges <- sorted_graph$g %>% igraph::get.edgelist()
  edges <- edges[nrow(edges):1, ]
  for (i in 1:nrow(edges)) {
    start <- edges[i, 1]
    end <- edges[i, 2]
    if (data_clus_cent[, 16][start] > data_clus_cent[, 16][end]) {
      edges[i, 1] <- end
      edges[i, 2] <- start
    }
  }
  edge_length <- igraph::E(sorted_graph$g)$weight %>% rev()

  milestone_network <- tibble::tibble(
    from = as.character(edges[, 1]),
    to = as.character(edges[, 2]),
    length = edge_length,
    directed = TRUE
  )

  out <- list(
    pseudotime = pseudotime,
    cluster = cl,
    data_clus_cent = data_clus_cent,
    milestone_network = milestone_network,
    g = sorted_graph
  )
  return(out)
}

#' @title goolam
#'
#' @description goolam datas set save as a SCE object.
"goolam"

#' @title genesets
#'
#' @description A list consists of pathway genesets of Homo sapiens (human) and Mus musculus (house mouse).
"genesets"

