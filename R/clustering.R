#' Affinity propagation clustering
#' 
#' AP clustering on genes and calculate the eigengenes of each cluster
#' @importFrom apcluster apcluster
#' @importFrom stats cor prcomp
#' @importFrom parallel makeCluster stopCluster
#' @importFrom S4Vectors metadata
#' @importFrom Matrix Matrix spMatrix colMeans t
#' @importFrom qlcMatrix corSparse cosSparse
#' @importFrom irlba prcomp_irlba
#' @import doParallel
#' @import foreach
#' @param sce A SingleCellExperiment object or matrix.
#' @param datatype Specify the data type in sce for filtering.
#' @param genes.use Specify the genes to perform clustering.
#' @param projections Whether to calculate reconstructions and scores using the L1 norm ("l1") the L2 norm ("l2").
#' @param sim Prior defined similarty matrix of genes for clustering - default is NULL.
#' @param method The method to calculate the affinity matrix, can be pearson (Defult) or cosine.
#' @param center Center the expression before performing PCA - default is FALSE.
#' @param scale. Scale the expression before performing PCA - default is FALSE.
#' @param normalize.score Wheter to normalize the pc1 score.
#' @param minsize Minimum cluster size.
#' @param ncore Number of cores used for parallel. 
#' @param ... Additional arguments passed on to \code{\link[apcluster]{apcluster}}
#' @return 
#' \itemize{
#'   \item Cluster labels of genes
#'   \item Basic stats of clusters
#'   \item Eigengenes
#' }
#' @export
ap.cluster <- function(sce, 
                       datatype = "logcounts",
                       genes.use = "genes.use",
                       projections = "l2",
                       sim = NULL,
                       method = "pearson",
                       center = TRUE, 
                       scale. = TRUE,
                       normalize.score = FALSE,
                       minsize = 2,
                       ncore = 2,
                       ...){
  # AP clustering genes and calculate the PCA pattern
  if (class(sce) != "SingleCellExperiment"){
    stop("sce must be a sce object")
  }
  if (!is.null(datatype)){
    if (!datatype %in% names(assays(sce))){
      stop("Available datatype are", names(assays(sce)))
    }
  }
  if (!is.null(genes.use)){
    if (!genes.use %in% colnames(rowData(sce))){
      stop("genes.use is not in rowData")
    }else{
      genes.use <- rowData(sce)[, genes.use]
    }
  }else{
    genes.use <- rep(TRUE, nrow(sce))
  }
  mat <- assay(sce[genes.use, ], datatype)
  if (is.null(sim)){
    if (method == "pearson"){
      sim <- as.matrix(corSparse(t(mat)))
    }
    else if (method == "cosine"){
      sim <- as.matrix(cosSparse(t(mat)))
    }
  }else{
    sim <- as.matrix(sim)
  }
  if (!"symbol" %in% colnames(rowData(sce))){
    rowData(sce) <- cbind(data.frame(symbol=rownames(sce)), as.data.frame(rowData(sce)))
  }
  dimnames(sim) <- list(rownames(mat), rownames(mat))
  message("Perform AP clustering on ", nrow(mat), " genes...")
  apres <- suppressWarnings(apcluster(sim, ...))
  labels <- lapply(apres@clusters, function(cl) data.frame(symbol = names(cl)))
  gene_num = sapply(labels,nrow)
  idx <- rep(1:length(labels), gene_num)
  labels <- do.call(rbind, labels)
  labels$labels <- idx
  labels$exemplars = as.vector(labels$symbol) %in% names(unlist(apres@exemplars))
  cls <- split(labels, labels$labels)
  cls[names(cls)[gene_num < minsize]] <- NULL
  used_cls = as.numeric(names(cls))
  message("Find ", length(used_cls), " clusters")

  ## L2 PCA ##
  if (projections == "l2"){
    if (ncore == 1){
      message("Calculate eigengenes...")
      pc1s <- do.call(cbind, lapply(cls, function(cl) prcomp_irlba(t(mat[as.character(cl$symbol),]), n=1, center = center, scale. = scale.)$x[,1]))
    }else{
      message("Calculate eigengenes using ", ncore, " cores...")
      cl <- makeCluster(ncore)
      registerDoParallel(cl, nocompile=FALSE)
      pc1s <- foreach(i = 1:length(cls), .packages=c("irlba"), .combine=cbind) %dopar% {
        prcomp_irlba(t(mat[as.character(cls[[i]]$symbol),]), n=1, center = center, scale. = scale.)$x[,1]
      }
      stopCluster(cl)
    }
  }
  ## L1 PCA ##
  else if (projections == "l1"){
    x <- "pcaL1"
    if (!require(x,character.only = TRUE)){
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop(paste0("Cannot install ", x))
    }
    if (ncore == 1){
      message("Perform wPCA and calculate eigengenes...")
      pc1s <- do.call(cbind, lapply(cls, function(cl) as.vector(pcaL1::wl1pca(t(as.matrix(mat[as.character(cl$symbol),])), projDim=1, center=center, projections=projections)$scores)))
    }
    else{
      message("Perform wPCA using ", ncore, " cores...")
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      pc1s <- foreach(i = 1:length(cls), .packages=c("pcaL1"), .combine=cbind) %dopar% {
        as.vector(wl1pca(t(as.matrix(mat[as.character(cls[[i]]$symbol),])), projDim=1, center=center, projections=projections)$scores)
      }
      stopCluster(cl)
    }
  }else{
    stop("projections should be l1 or l2")
  }
  m <- do.call(cbind,lapply(cls, function(cl) colMeans(mat[as.character(cl$symbol),])))

  # Correct the expression direction of PC1
  direction <- diag(sign(cor(pc1s, y = m)))
  pc1s <- sweep(pc1s, 2, direction, "*")
  colnames(pc1s) <- paste0("cl", names(cls)) #
  rownames(pc1s) <- colnames(mat)
  ap_stats <- data.frame(id = paste0("cl",used_cls),
                         cluster_gene = paste0("cl",used_cls,"_",names(unlist(apres@exemplars))[used_cls]),
                         cluster_exemplar = names(unlist(apres@exemplars))[used_cls],
                         cluster_size = gene_num[used_cls])
  
  # Plot cluster size
  df <- ap_stats
  df$color <- as.factor(rep(c(1,2),(nrow(df)+1)/2)[1:nrow(df)])
  p <- ggplot(data=df, aes(x=factor(cluster_gene, levels = df$cluster_gene), y=cluster_size, fill=color)) +
        geom_bar(stat="identity",width = 0.8) + scale_fill_manual(values=c("#3DB9BF","#E87167"))+
        labs(title = "Cluster size", x = "Cluster exemplars", y = "Gene number") +
        theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size=6), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none") 
  metadata(sce)$p.cluster.size <- p
  message("Use metadata(sce)$p.cluster.size to visualize cluster size.")
  
  labels <- labels[labels$labels %in% used_cls, ]
  lv1.ap <- labels[match(rownames(sce), labels$symbol),]
  rowData(sce)[, c("lv1.labels", "lv1.exemplars")] <- lv1.ap[, c("labels", "exemplars")]
  if (normalize.score){
    reducedDim(sce, "lv1.score") <- scale(pc1s)
  }else{
    reducedDim(sce, "lv1.score") <- pc1s
  }
  metadata(sce)$lv1.parameters <- c("projections"=projections, "method"=method)
  metadata(sce)$ap.stats <- ap_stats
  return(sce)
}

#' Reduce redundant clusters of AP clustering 
#' 
#' Merge highly related clusters into groups
#' @importFrom apcluster apcluster
#' @importFrom stats cor prcomp
#' @importFrom parallel makeCluster stopCluster
#' @importFrom S4Vectors metadata
#' @importFrom qlcMatrix corSparse cosSparse
#' @importFrom Matrix Matrix spMatrix colMeans t
#' @importFrom irlba prcomp_irlba
#' @import doParallel
#' @import foreach
#' @param sce A SingleCellExperiment object.
#' @param datatype Sepcify the data type in sce for filtering.
#' @param method The method to calculate the affinity matrix, can be pearson (Defult) or cosine.
#' @param dsim Prior defined similarty matrix of lv1.score for clustering - default is NULL.
#' @param center Center the expression before performing PCA - default is FALSE.
#' @param scale. Scale the expression before performing PCA - default is FALSE.
#' @param normalize.score Whether to normalize pc1 score.
#' @param ncore Number of cores used for parallel. 
#' @param ... Additional arguments passed on to \code{\link[apcluster]{apcluster}}
#' @return 
#' \itemize{
#'   \item Group labels of genes
#'   \item Basic stats of Groups
#'   \item Eigengenes
#' }
#' @export
reduce.cluster <- function(sce, 
                           datatype = "logcounts",
                           genes.use = "genes.use",
                           projections = "l2",
                           method = "pearson",
                           sim = NULL,
                           center = TRUE, 
                           scale. = TRUE,
                           normalize.score = TRUE,
                           ncore = 2,
                           ...){

  if (class(sce) != "SingleCellExperiment"){
    stop("sce must be a sce object")
  }
  if (!is.null(datatype)){
    if (!datatype %in% names(assays(sce))){
      stop("Available datatype are", names(assays(sce)))
    }
  }
  if (!is.null(genes.use)){
    if (!genes.use %in% colnames(rowData(sce))){
      stop("genes.use is not in rowData")
    }else{
      genes.use <- rowData(sce)[, genes.use]
    }
  }else{
    genes.use <- rep(TRUE, nrow(sce))
  }
  
  mat = assay(sce[genes.use, ], datatype)
  if (!"symbol" %in% colnames(rowData(sce))){
    rowData(sce) <- cbind(data.frame(symbol=rownames(sce)), as.data.frame(rowData(sce)))
  }
  gc <- rowData(sce)[,c("symbol","lv1.labels","lv1.exemplars")]
  gc <- gc[complete.cases(gc), ]
  ap_stats <- metadata(sce)$ap.stats
  y <- reducedDim(sce,"lv1.score")
  if (is.null(sim)){
    if (method == "pearson"){
      sim <- as.matrix(corSparse(y))
    }
    else if (method == "cosine"){
      sim <- as.matrix(cosSparse(y))
    }
  }else{
    sim <- as.matrix(sim)
  }
  dimnames(sim) <- list(colnames(y), colnames(y))
  message("Perform AP clustering on ", ncol(y), " clusters...")
  apres <- suppressWarnings(apcluster(sim, ...))
  r2 <- lapply(apres@clusters, function(cl) data.frame(cluster = names(cl)))

  idx <- rep(1:length(r2), sapply(r2,nrow))
  r2 <- do.call(rbind, r2)
  r2$groups <- idx

  r2l <- split(r2$cluster, r2$groups)
  grl <- lapply(r2l, function(gr){
    al <- gc[gc$lv1.labels %in% as.integer(gsub("*.([0-9]*)","\\1",gr)),]
    al[!duplicated(as.character(al$symbol)), ]
  })

  gr_num <- sapply(grl, nrow)
  idx <- rep(1:length(r2l), gr_num)
  ndf <-  Reduce(rbind, grl)
  ndf$group <- idx # Add group labels

  ## L2 PCA ##
  if (projections == "l2"){
    if (ncore == 1){
      message("Calculate eigengenes...")
      pc1s <- do.call(cbind, lapply(grl, function(gr) prcomp_irlba(t(mat[as.character(gr$symbol),]), n = 1, center = center, scale. = scale.)$x[,1]))
    }else{
      message("Calculate eigengenes using ", ncore, " cores...")
      cl <- makeCluster(ncore)
      registerDoParallel(cl, nocompile=FALSE)
      pc1s <- foreach(i = 1:length(grl), .packages=c("irlba"), .combine=cbind) %dopar% {
        prcomp_irlba(t(mat[as.character(grl[[i]]$symbol),]), n = 1, center = center, scale. = scale.)$x[,1]
      }
      stopCluster(cl)
    }

  }
  ## L1 PCA ##
  else if (projections == "l1"){
    x <- "pcaL1"
    if (!require(x,character.only = TRUE)){
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop(paste0("Cannot install ", x))
    }
    if (ncore == 1){
      message("Perform wPCA and calculate eigengenes...")
      pc1s <- do.call(cbind, lapply(grl, function(gr) as.vector(pcaL1::wl1pca(t(as.matrix(mat[as.character(gr$symbol),])), projDim=1, center=FALSE, projections=projections)$scores)))
    }
    else{
      message("Perform wPCA using ", ncore, " cores...")
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      pc1s <- foreach(i = 1:length(grl), .packages=c("pcaL1"), .combine=cbind) %dopar% {
        as.vector(wl1pca(t(as.matrix(amat[as.character(grl[[i]]$symbol),])), projDim=1, center=FALSE, projections=projections)$scores)
      }
      stopCluster(cl)
    }
  }
  message("Find ", length(grl), " groups")
  
  # Correct the expression direction of PC1
  m <-  do.call(cbind,lapply(grl, function(gr) colMeans(mat[as.character(gr$symbol),])))
  direction <- diag(sign(cor(pc1s, y = m)))
  pc1s <- sweep(pc1s, 2, direction, "*")
  colnames(pc1s) <- paste0("group", names(grl))
  rownames(pc1s) <- colnames(mat)
  r2 <- r2[match(ap_stats$id, r2$cluster),]
  ap_stats$group <- r2$groups
  ap_stats$group_size <- gr_num[ap_stats$group]
  
  # Plot group size
  df <- ap_stats[order(ap_stats$group),]
  df <- df[!duplicated(df$group),]
  df$color <- as.factor(rep(c(1,2),(nrow(df)+1)/2)[1:nrow(df)])
  p <- ggplot(data=df, 
              aes(x=factor(group, levels = df$group), y=group_size, fill=color)) +
    geom_bar(stat="identity",width = 0.8) + scale_fill_manual(values=c("#3DB9BF","#E87167"))+
    labs(title = "Group size", x = "Group", y = "Gene number") +
    theme_bw() + theme(axis.text.x = element_text(size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none") 
  metadata(sce)$p.group.size <- p
  message("Use metadata(sce)$p.group.size to visualize group size.")
  
  if (normalize.score){
    reducedDim(sce, "lv2.score") <- scale(pc1s)
  }else{
    reducedDim(sce, "lv2.score") <- pc1s
  }
  lv2.ap <- ndf[match(rownames(sce), ndf$symbol),]
  rowData(sce)$lv2.labels <- lv2.ap[,"group"]
  metadata(sce)$lv2.parameters <- c("projections"=projections, "method"=method)
  metadata(sce)$ap.stats <- ap_stats
  
  return(sce)
  
}

#' Heatmap for gene expression pattern.
#' 
#' Plot the clustering result of cells and genes.
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @importFrom utils head
#' @param mat Expression matrix, columns are cells and rows are genes.
#' @param filename Filename to save your heatmap.
#' @param vmax The maximum value to plot.
#' @param vmin The minimum value to plot.
#' @param log log2 tranform the data.
#' @param scale Scale the data, result in mean = 0 and std = 1 of each gene.
#' @param hlight The genes you want to highlight in the heatmap.
#' @param cmap The pallette of the heatmap, one can generate it using colorRampPalette().
#' @param cluster_rows boolean values determining if rows should be clustered or hclust object.
#' @param cluster_cols boolean values determining if columns should be clustered or hclust object.
#' @param gene_labels A vector of gene clustering labels.
#' @param cell_labels A vector of cell clustering labels.
#' @param gene_order Specify a sequence of gene clusters order, the values must be in gene_labels.
#' @param cell_order Specify a sequence of cell clusters order, the values must be in cell_labels.
#' @param gaps_row A vector indicates gaps of row, if NULL the gaps will be setting based on the gene_labels.
#' @param gaps_col A vector indicates gaps of column, if NULL the gaps will be setting based on the cell_labels.
#' @param anno_row Data frame that specifies the annotations shown on left side of the heatmap. Each row defines the features for a specific row. The rows in the data and in the annotation are matched using corresponding row names. Note that color schemes takes into account if variable is continuous or discrete.
#' @param anno_col Similar to annotation_row, but for columns.
#' @param plot_gaps_row Make gaps on row.
#' @param plot_gaps_col Make gaps on column.
#' @param re_write_gene_labels Re-write gene labels to an ordered sequeence.
#' @param re_write_cell_labels Re-write cell labels to an ordered sequeence.
#' @param ... Additional arguments passed on to pheatmap.
#' @return Ordering cluster results.
#' @export
cheatmap <- function(mat,
                     filename = NA,
                     vmax = 2,
                     vmin = -1, 
                     log = FALSE,
                     scale = FALSE,
                     hlight = NULL,
                     cmap = NULL,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     gene_labels = NULL,
                     cell_labels = NULL,
                     gene_order =NULL,
                     cell_order =NULL, 
                     gaps_row = NULL,
                     gaps_col = NULL,
                     anno_row = NULL,
                     anno_col = NULL,
                     plot_gaps_row = TRUE,
                     plot_gaps_col = TRUE,
                     re_write_gene_labels = FALSE,
                     re_write_cell_labels = FALSE,
                     show_colnames = FALSE,
                     ...){
  
  if (!is.numeric(vmax) || !is.numeric(vmin)){
    stop("The expression boundary must be numeric")
  }
  if (length(gene_order) == 1){
    gene_order <- as.numeric(strsplit(as.character(gene_order),"")[[1]])
  }
  if (length(cell_order) == 1){
    cell_order <- as.numeric(strsplit(as.character(cell_order),"")[[1]])
  }
  if (!is.null(gene_order) && !is.null(gene_labels) && length(intersect(gene_order, unique(gene_labels)))!=length(unique(gene_labels))){
    stop("The elements of gene_order must be consistent with gene_labels: ",unique(gene_labels))
  }
  if (!is.null(cell_order) && !is.null(cell_labels) && length(intersect(cell_order, unique(cell_labels)))!=length(unique(cell_labels))){
    stop("The elements of cell_order must be consistent with cell_labels: ",unique(cell_labels))
  }
  
  # Order genes
  makeGeneOrder <- 1:nrow(mat)
  if (!is.null(gene_labels)){
    cluster_rows <- FALSE
    if (!is.null(gene_order)){
      makeGeneOrder <- alter_label(gene_labels, gene_order)
      if (plot_gaps_row){
        if (is.null(gaps_row)){
          gaps_row <- head(cumsum(table(gene_labels)[as.character(gene_order)]),-1)
        }
      }
      if (re_write_gene_labels){
        gene_labels <- re_order(gene_labels, gene_order)
      }
    }
    else {
      makeGeneOrder <- order(gene_labels)
      if (plot_gaps_row){
        gaps_row <- head(cumsum(table(gene_labels)),-1)
      }
    }
    mat <- mat[makeGeneOrder,]
    output_gene_labels <- data.frame(labels=gene_labels[makeGeneOrder],row.names = rownames(mat))
    if (is.null(anno_row)){
      anno_row = data.frame(Module=factor(output_gene_labels$labels, levels = unique(output_gene_labels$labels)),
                            row.names = rownames(mat))
    }
  }
  else{
    gaps_row <- NULL
    output_gene_labels <- NULL
  }
  
  # Order cells
  makeCellOrder <- 1:ncol(mat)
  if (!is.null(cell_labels)){
    cluster_cols <- FALSE
    if (!is.null(cell_order)){
      makeCellOrder <- alter_label(cell_labels, cell_order)
      if (plot_gaps_col){
        if (is.null(gaps_col)){
          gaps_col <- head(cumsum(table(cell_labels)[as.character(cell_order)]),-1)
        }
      }
      if (re_write_cell_labels){
        cell_labels <- re_order(cell_labels, cell_order)
      }
    }
    else {
      makeCellOrder <- order(cell_labels)
      if (plot_gaps_col){
        gaps_col <- head(cumsum(table(cell_labels)),-1)
      }
    }
    mat <- mat[,makeCellOrder]
    output_cell_labels <- data.frame(labels=cell_labels[makeCellOrder],row.names = colnames(mat))
    if (is.null(anno_col)){
      anno_col <- data.frame(Cluster = factor(output_cell_labels$labels, levels = unique(output_cell_labels$labels)), 
                             row.names = colnames(mat)) # bestcells = bestcells EGFP=EGFP
    }
  }
  else{
    gaps_col <- NULL
    output_cell_labels <- NULL
  }
  
  if (log){
    mat <- log2(mat + 1)
  }
  if (scale){
    mat <- t(scale(t(mat)))
  }
  
  mat[mat>vmax] <- vmax
  mat[mat<vmin] <- vmin
  
  # hightlight genes
  psuGene = rep(c(""),dim(mat)[1])
  
  # hlight = tfsinRep
  if (!is.null(hlight)){
    for (g in hlight){
      idx = which(rownames(mat) == g)
      psuGene[idx] = g
    }
  }
  
  if (is.null(cmap)){
    cmap <- colorRampPalette(c("#7431A6","#000000","#F2F208"), bias = 1.7)(255)
  }
  pheatmap(mat,
           filename = filename,
           annotation_row = anno_row,
           annotation_col = anno_col,
           gaps_col = gaps_col,
           gaps_row = gaps_row,
           color = cmap,
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           labels_row = psuGene,
           show_colnames = show_colnames,
           ...)
  
  return(list(heat=mat,cell_labels=output_cell_labels,gene_labels=output_gene_labels))
  
}

#' Re-write clustering labels.
#' 
#' Re-write clustering labels.
#' @param labs Clustering labels.
#' @param ord_id Cluster order idx.
#' @return A vector of re-write labels.
#' @export
#' 
re_order <- function(labs,ord_id){
  new_ord <- 1:max(labs)
  names(new_ord) <- ord_id
  new_labs <- new_ord[as.character(labs)]
  
  return(new_labs)
  
}

#' Re-ordering clusters.
#' 
#' Re-ordering clusters for heatmap.
#' @param ori_labels Origin labels.
#' @param re_labels Re-ordered labels.
#' @return New labels after ordering.
#' @export
alter_label <- function(ori_labels,re_labels){
  if (length(intersect(unique(ori_labels),unique(re_labels)))!=length(unique(ori_labels))){
    stop("The given labels are not consistent with the origin labels")
  }
  in_order=seq(length(ori_labels))
  names(in_order) <- ori_labels
  tre_labelspe=sort(unique(ori_labels))
  l <- list()
  for(i in 1:length(tre_labelspe)){
    l[[i]]=in_order[which(names(in_order) %in% tre_labelspe[i])]
  }
  names(l) = tre_labelspe
  new_order=c()
  for(i in re_labels){
    new_order <- c(new_order,l[[as.character(i)]])
  }
  
  return(new_order)
  
}
