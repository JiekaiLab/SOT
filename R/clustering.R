#' Affinity propagation clustering
#' 
#' AP clustering on genes and calculate the eigengenes of each cluster
#' @importFrom apcluster apcluster
#' @importFrom stats cor prcomp
#' @param sce A SingleCellExperiment object or matrix.
#' @param datatype Sepcify the data type in sce for filtering.
#' @param method The method to calculate the affinity matrix.
#' @param power Power level of the paire-wised correlation.
#' @param aff Prior defined affinity matrix - default is NULL.
#' @param trim Trim the largest and lowest values of each gene to avoid outlier effect.
#' @param center Center the expression before performing PCA - default is FALSE.
#' @param scale. Scale the expression before performing PCA - default is FALSE.
#' @param ... Additional arguments passed on to \code{\link[apcluster]{apcluster}}
#' @return 
#' \itemize{
#'   \item Cluster labels of genes
#'   \item Basic stats of clusters
#'   \item Eigengenes
#' }
#' @export
ap.cluster <- function(sce, 
                       datatype = NULL,
                       method = "cor",
                       power = 6,
                       aff = NULL,
                       trim = 5.1/ncol(sce),  
                       center = TRUE, 
                       scale. = TRUE,
                       minsize = 5,
                       ...){
  # AP clustering genes and calculate the PCA pattern
  if (class(sce) == "SingleCellExperiment"){
    mat = assay(sce, datatype)
  }
  else(
    mat = as.matrix(sce)
  )
  if (trim > 0){
    mat <- scde::winsorize.matrix(mat, trim = trim) # Remove outliers with heighest and lowest expression
    vi <- which(abs(apply(mat, 1, function(x) sum(abs(diff(x))))) > 0)
    mat <- mat[vi, ]
  }
  
  if (is.null(aff)){
    if (method == "cor"){
      aff <- cor(t(mat))
    }
    else if(method == "wgcna"){
      message("Calculate topology overlap matrix")
      if ("WGCNA" %in% rownames(installed.packages())){
        aff <- WGCNA::TOMsimilarity(WGCNA::adjacency(t(mat), type = "signed", power = power))
      }
      else{
        aff <- TOMsimilarity(adjacency(t(mat), type = "signed", power = power))
      }
    }
    else{
      stop("method should be 'cor' or 'wgcna'")
    }
  }
  dimnames(aff) <- list(rownames(mat), rownames(mat))
  apres <- apcluster(aff, ...)
  cls <- apres@clusters
  gene_num = sapply(cls,length)
  idx <- rep(1:length(apres@exemplars), gene_num)
  labels <- data.frame(labels = idx, 
                       exemplars = unlist(apres@clusters) %in% unlist(apres@exemplars),
                       row.names = names(unlist(apres@clusters)))
  pca <- lapply(cls, function(cl) prcomp(t(as.matrix(mat[cl,])), center = center, scale. = scale.))
  var <- sapply(pca, function(cl) cl$sdev[1]^2)
  pc1s <- do.call(rbind, lapply(pca, function(cl) cl$x[,1]))
  # Correct the expression direction of PC1
  m <-  do.call(rbind,lapply(cls, function(cl) colMeans(mat[cl,])))
  direction <- diag(sign(cor(t(pc1s), y = t(m))))
  pc1s <- sweep(pc1s, 1, direction, "*")
  rownames(pc1s) <- paste0("cl", 1:nrow(pc1s))
  cl_stats <- data.frame(exemplars = names(unlist(apres@exemplars)),
                         var = var, 
                         gene_num = gene_num,
                         row.names = rownames(pc1s))
  
  cl_stats <- cl_stats[cl_stats$gene_num >= minsize, ]
  pc1s <- pc1s[rownames(cl_stats), ]
  labels <- labels[labels$labels %in% as.integer(gsub("*.([0-9]*)","\\1",rownames(cl_stats))), ]
  
  return(list(cluster = labels, cl_stats = cl_stats, y = pc1s))
  
}

#' Reduce redundant clusters of AP clustering 
#' 
#' Merge highly related clusters into groups
#' @importFrom apcluster apcluster
#' @importFrom stats cor prcomp
#' @param sce A SingleCellExperiment object or matrix.
#' @param datatype Sepcify the data type in sce for filtering.
#' @param gcl AP clustering result for the firt round.
#' @param aff Prior defined affinity matrix - default is NULL.
#' @param trim Trim the largest and lowest values of each gene to avoid outlier effect.
#' @param center Center the expression before performing PCA - default is FALSE.
#' @param scale. Scale the expression before performing PCA - default is FALSE.
#' @param ... Additional arguments passed on to \code{\link[apcluster]{apcluster}}
#' @return 
#' \itemize{
#'   \item Group labels of genes
#'   \item Basic stats of Groups
#'   \item Eigengenes
#' }
#' @export
reduce.cluster <- function(sce, 
                           gcl,
                           datatype = NULL,
                           aff = NULL,
                           trim = 5.1/ncol(mat),
                           center = TRUE, 
                           scale. = TRUE,
                           ...){
  if (class(sce) == "SingleCellExperiment"){
    mat = assay(sce, datatype)
  }
  else(
    mat = as.matrix(sce)
  )
  gc <- gcl$cluster
  gs <- gcl$cl_stats
  y <- gcl$y
  
  # AP group pc1 pattern of gene clusters
  r2 <- ap.cluster(y, aff = aff, trim = trim, ...)
  r2l <- split(r2$cluster, r2$cluster$labels)
  grl <- lapply(r2l, function(gr){
    gc[gc$labels %in% as.integer(gsub("*.([0-9]*)","\\1",rownames(gr))),]
  })
  
  gr_num <- sapply(grl, nrow)
  idx <- rep(1:length(r2l), gr_num)
  ndf <-  Reduce(rbind, grl)
  ndf$group <- idx # Add group labels
  
  # Calculate PC1 score of groups
  if (trim > 0){
    mat <- scde::winsorize.matrix(as.matrix(mat), trim = trim) # Remove outliers with heighest and lowest expression
    vi <- which(abs(apply(mat, 1, function(x) sum(abs(diff(x))))) > 0)
    mat <- mat[vi, ]
  }
  pca <- lapply(grl, function(gr) prcomp(t(as.matrix(mat[rownames(gr),])), center=center, scale. = scale.))
  var <- sapply(pca, function(gr) gr$sdev[1]^2)
  pc1s <- do.call(rbind, lapply(pca, function(gr) gr$x[,1]))
  
  # Correct the expression direction of PC1
  m <-  do.call(rbind,lapply(grl, function(gr) colMeans(mat[rownames(gr),])))
  direction <- diag(sign(cor(t(pc1s), y = t(m))))
  pc1s <- sweep(pc1s, 1, direction, "*")
  rownames(pc1s) <- paste0("group", 1:nrow(pc1s))
  gs$group <- r2$cluster[rownames(gs),]$labels
  gs$group_var <- var[gs$group]
  gs$group_num <- gr_num[gs$group]
  
  return(list(group = ndf, gr_stats = gs, y = pc1s))
  
}
##
