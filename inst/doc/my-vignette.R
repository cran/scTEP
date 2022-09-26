## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("scTEP")

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages({
library(SummarizedExperiment)
library(scTEP)
})

## ---- eval=FALSE--------------------------------------------------------------
#  #Load example data (SCE dataset)
#  data("goolam")
#  #Get data matrix and label
#  expr <- as.matrix(t(assay(goolam)))
#  label <- as.character(goolam$label)
#  stages = goolam@metadata$cell.stages

## ---- eval=FALSE--------------------------------------------------------------
#  dim(expr)
#  expr[1:10,1:10]

## ---- eval=FALSE--------------------------------------------------------------
#  data = preprocessing(expr)
#  data$expr[1:10,1:10]

## ---- eval=FALSE--------------------------------------------------------------
#  data("genesets")
#  genesets$mmu$`path:mmu00010`
#  
#  data_fa = scTEP.fa(data, genesets, data_org = 'mmu', seed = 1)
#  dim(data_fa)
#  data_fa[1:10,1:10]

## ---- eval=FALSE--------------------------------------------------------------
#  allCluster = scTEP::clustering(data, seed = 1)

## ---- eval=FALSE--------------------------------------------------------------
#  scDHA_res <- scDHA(data_fa, do.clus = T, gen_fil = T, ncores = 16, seed = 1)

## ---- eval=FALSE--------------------------------------------------------------
#  idx = which(label == stages[1])
#  out = trajectoryinference(data, idx, scDHA_res, allCluster, seed = 1)
#  plot(out$g$g)
#  r = round(cor(out$pseudotime, as.numeric(factor(label, levels = stages))), digits = 2)

## ---- eval=FALSE--------------------------------------------------------------
#  suppressPackageStartupMessages({
#  library(irlba)
#  library(ggplot2)
#  library(uwot)
#  })
#  # Get 2D emebedding data of the original data.
#  cols <- c(
#      "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF",
#      "#EE4C97FF", "#000075", "#a9a9a9", "#8DD3C7", "#C8EABC", "#FBFBB4", "#D9D7C9", "#C3B4D0",
#      "#E39699", "#E9877F", "#A9A0B2", "#97B1BD", "#D9B382",
#      "#EBBD63", "#C4D367", "#C7D98C", "#EED0CD", "#F0D1E1",
#      "#DED7DA", "#CDB7CE", "#BE88BF", "#C2ADC0", "#CBE5C4",
#      "#E4EB9C", "#FFED6F", "#CCFF99", "#33FF00", "#FFDB6D", "#33CC33", "#003300", "#00CC99",
#      "#FFFF00", "#CC9900", "#FFFFCC", "#CCFFCC",
#      "#2059BB", "#16489E", "#0F3980", "#0E2F68", "#8A3ABF",
#      "#7330A0", "#5A3870", "#DE1A64", "#C21A59", "#6D1234", "#3D3135", "#2ABAA4", "#5C9F95",
#      "#335650", "#D55C31", "#B07966", "#E2D8D4"
#    )
#  umap.res <- uwot::umap(latent)

## ---- eval=FALSE--------------------------------------------------------------
#  set.seed(1)
#  cell_stages <- factor(label, levels = stages)
#  label <- as.numeric(factor(label, levels = stages))
#  latent <- cbind(umap.res, label) %>% as.data.frame()
#  
#  # Calculate cluster centroid and create dataframe saving start and end points.
#  clus_center <- lapply(1:length(unique(label)), function(cl) colMeans(latent[label == cl, ])) %>%
#    do.call(what = rbind)
#  colnames(clus_center) <- c("x", "y", "cluster_id")
#  
#  figure <- ggplot2::ggplot(latent, ggplot2::aes(x = V1, y = V2, color = cell_stages)) +
#    ggplot2::geom_point() +
#    # Plot the centroid
#    ggplot2::geom_point(data = as.data.frame(clus_center), ggplot2::aes(x = x, y = y), color = 'black', size = 2) +
#    ggplot2::labs(x = paste0("UMAP1"), y = paste0("UMAP2"), title = paste("Landscape")) +
#    ggplot2::theme_classic() +
#    ggplot2::scale_color_manual(values = cols) +
#    # annotate(geom = "table", x=min(tsne_original$V1),y=min(tsne_original$V2), label = silhou)+
#    # scale_color_npg()+
#    ggsci::scale_fill_npg() +
#    ggplot2::theme(
#      legend.position = "right"
#    )
#  figure

## ---- eval=FALSE--------------------------------------------------------------
#  latent <- uwot::umap(scDHA_res$latent) %>% as.data.frame()
#  pseudotime = out$pseudotime
#  # Calculate cluster centroid and create dataframe saving start and end points.
#  figure <- ggplot2::ggplot(latent, ggplot2::aes(x = V1, y = V2)) +
#    ggplot2::labs(x = paste0("UMAP1"), y = paste0("UMAP2"), title = paste("Pseudotime landscape")) +
#    ggplot2::geom_point(ggplot2::aes(color = pseudotime), size = 2) +
#    ggplot2::theme_classic() +
#    ggplot2::theme(
#      legend.position = "right",
#      plot.title = ggplot2::element_text(),
#      axis.text = ggplot2::element_text(),
#      axis.title = ggplot2::element_text(),
#      legend.key.size = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(), legend.title = ggplot2::element_text)
#  figure

## ---- eval=FALSE--------------------------------------------------------------
#  # Draw trajectory plot
#  edge_nodes <- out$milestone_network[, 1:2] %>% as.data.frame()
#  cluster = out$cluster
#  latent <- cbind(umap.res, cluster) %>% as.data.frame()
#  cell_stages <- factor(cluster)
#  
#  # Calculate cluster centroid and create dataframe saving start and end points.
#  clus_center <- lapply(1:length(unique(cluster)), function(cl) colMeans(latent[which(cluster == cl), ])) %>%
#    do.call(what = rbind)
#  
#  colnames(clus_center) <- c("x", "y", "cluster_id")
#  clulines <- NULL
#  for (edge_index in 1:nrow(edge_nodes)) {
#    start_cluster <- edge_nodes[edge_index, 1] %>% as.numeric()
#    end_cluster <- edge_nodes[edge_index, 2] %>% as.numeric()
#    temp_points <- c(
#      clus_center[, "x"][start_cluster], clus_center[, "y"][start_cluster],
#      clus_center[, "x"][end_cluster], clus_center[, "y"][end_cluster]
#    )
#    clulines <- rbind(clulines, temp_points)
#  }
#  colnames(clulines) <- c("x", "y", "xend", "yend")
#  clulines <- as.data.frame(clulines)
#  clulines$x.mid <- (clulines$x + clulines$xend) / 2
#  clulines$y.mid <- (clulines$y + clulines$yend) / 2
#  
#  
#  figure <- ggplot2::ggplot(latent, ggplot2::aes(x = V1, y = V2, color = cell_stages)) +
#    ggplot2::geom_point() +
#    # Plot the centroid
#    ggplot2::geom_point(data = as.data.frame(clus_center), ggplot2::aes(x = x, y = y), color = "black", size = 2) +
#    ggplot2::labs(x = paste0("UMAP1"), y = paste0("UMAP2"), title = paste("Trajectory")) +
#    # Add Sting to the figure
#    ggplot2::geom_segment(ggplot2::aes_string(x = "x.mid", xend = "xend", y = "y.mid", yend = "yend", size = NULL),
#                          data = clulines, color = "black", size = 1.25
#    ) +
#    # Add arrow to string
#    ggplot2::geom_segment(
#      arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"), type = "closed", ends = "last"),
#      ggplot2::aes_string(x = "x", xend = "x.mid", y = "y", yend = "y.mid", size = NULL),
#      data = clulines, color = "black", size = 1.25
#    ) +
#    ggplot2::theme_classic() +
#    ggplot2::scale_color_manual(values = cols) +
#    ggsci::scale_fill_npg() +
#    ggplot2::theme(
#      legend.position = "right",
#      plot.title = ggplot2::element_text(),
#      axis.text = ggplot2::element_text(),
#      axis.title = ggplot2::element_text(),
#      legend.key.size = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(), legend.title = ggplot2::element_text()
#    )
#  figure

## ---- eval=FALSE--------------------------------------------------------------
#  #Draw_dev_pseudotime <- function(dataset, label, cell.stages, pseudotime, r, cols) {
#  raw <- goolam$label %>% as.character() %>% as.data.frame()
#  raw$time <- pseudotime
#  raw[, 1] <- factor(raw[, 1], levels = stages)
#  colnames(raw) <- c("cell_type","pseudotime")
#  figure <- ggplot2::ggplot(raw, ggplot2::aes(y = cell_type, x = pseudotime, color = cell_type)) +
#    ggplot2::geom_jitter() +
#    ggplot2::labs(x = "Pseudo Time", y = "", title = paste0("Pseudotime")) +
#    ggplot2::scale_color_manual(values = cols) +
#    ggplot2::theme_classic() +
#    ggplot2::theme(
#      legend.position = "right",
#      plot.title = ggplot2::element_text(size = 20),
#      axis.text = ggplot2::element_text(size = 20),
#      axis.title = ggplot2::element_text(size = 20),
#      legend.key.size = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size = 20), legend.title = ggplot2::element_text(size = 20),
#    )
#  figure

