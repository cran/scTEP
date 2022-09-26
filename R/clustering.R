#' @importFrom scDHA scDHA
#' @importFrom parallel makeCluster clusterEvalQ stopCluster
#' @importFrom foreach %dopar% foreach
#' @importFrom rlang .data
#' @title scTEP
#' @description The 'clustering' function conducts multiple clustering to the sc-RNA seq data using the 'scDHA' function from the 'scDHA' package. The 'scDHA' allows the user to set a specific cluster number for the clustering process. We set the cluster number from 6 to 10 and run the 'scDHA' function six times. The 'trajectoryinference' function will using those total six clustering results to generate pseudotimes.
#' @param data A list consists of gene expression matrix.
#' @param ncores Number of processor cores to use. This values is set to \code{seed = 10L} by default.
#' @param seed A parameter to set a seed for reproducibility.
#' @return List with the following keys:
#' \itemize{
#' \item allCluster - A list consists of clutering results using scDHA with \code{k = 5:10}.
#' }
#'
#' @references
#'
#' 1. Duc Tran, Hung Nguyen, Bang Tran, Carlo La Vecchia, Hung N. Luu, Tin Nguyen (2021). Fast and precise single-cell data analysis using a hierarchical autoencoder. Nature Communications, 12, 1029. doi: 10.1038/s41467-021-21312-2
#'
#' @examples
#' \donttest{
#' # Load the package and the example data (goolam data set)
#' library(scTEP)
#' #Load pathway genesets
#' data('genesets')
#' #Load example data (SCE dataset)
#' data("goolam")
#' #Get data matrix and label
#' expr <- as.matrix(t(SummarizedExperiment::assay(goolam)))[, 1:100]
#'
#' #Get data matrix and label
#' data = preprocessing(expr)
#'
#' #Get clustering results using scDHA with k from 6 to 10.
#' allCluster = scTEP::clustering(data, ncores = 2)
#' }
#' @export
clustering <- function(data, ncores = 10L, seed = NULL) {
  cl <- parallel::makeCluster(min(3, ncores), outfile = "/dev/null")
  doParallel::registerDoParallel(cl, cores = ncores)
  parallel::clusterEvalQ(cl, {
    library(tidyverse)
    library(scDHA)
    library(matrixStats)
    library(doParallel)
    library(SummarizedExperiment)
    library(SingleCellExperiment)
    library(Matrix)
    library(foreach)
  })
  k <- NULL
  allCluster <- foreach::foreach(k = 5:10) %dopar% {
    tmp <- (function() {
      scDHARes <- scDHA::scDHA(data$expr, k = k, do.clus = T, gen_fil = T, ncores = ncores, seed = seed)
      scDHARes$cluster
    })()
  }
  parallel::stopCluster(cl)
  allCluster <- lapply(allCluster, as.integer)
  return(allCluster)
}
