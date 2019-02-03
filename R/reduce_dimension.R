#' Force-directed layout
#' 
#' 2-D visualization of data on kNN graph using force-directed layout
#' @importFrom FNN knn.index
#' @importFrom igraph graph_from_adj_list layout_with_fr degree cluster_louvain cluster_spinglass as.undirected
#' @importFrom stats cor prcomp
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom SingleCellExperiment reducedDim reducedDimNames colData
#' @importFrom SummarizedExperiment assay colData
#' @importFrom utils installed.packages install.packages
#' @importFrom ggpubr ggarrange rremove
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @param sce A SingleCellExperiment object.
#' @param usedEmbed Specify the reduce dimension type for construct network. See current avaible embedding by \code{\link[SingleCellExperiment]{reducedDimNames}}.
#' @param usedDims A vector to specify particular dimensions in usedEmbed to calculate.
#' @param datatype Type of expression.
#' @param filename Filename to save - default is NULL.
#' @param k k nearest neighbors to construct graph.
#' @param n_trees Number of trees for largeVis.
#' @param max_iter max_iter for largeVis. see \code{\link[SOT]{randomProjectionTreeSearch}}
#' @param seed Random seed.
#' @param layout_method Graph layout method of FR (Fruchterman-Reingold layout) or FA2 (ForceAtlas2).
#' @param grid Character scalar, whether to use the faster, but less accurate grid based implementation of the algorithm. By default (“auto”), the grid-based implementation is used if the graph has more than one thousand vertices.
#' @param iterations Iterations number of iterations to be performed.
#' @param clus_method Clustering method besed on graph.
#' @param colour_condition Label the cells with continuous color.
#' @param repel The repel constant: the greater the constant k the stronger the repulsion force between points.
#' @param plotstep Plot step.
#' @param colour_gene A character of a vector of a gene expression.
#' @param color_scale Scatter plot color correponding to levels - default is NULL.
#' @param colour_embed A character of a reduced dimentsion.
#' @param cmap Color style of cmap will be used when color_scale is NULL.
#' @param compute_layout Compute fr layout or not. If compute_layout is FASLE, uese the usedEmbed as layout.
#' @param orderExp Order plot to highlight high expressed cells.
#' @param title Plot title.
#' @param ncol (optional) number of columns in the plot grid.
#' @param nrow (optional) number of rows in the plot grid.
#' @param alpha Transparency of points.
#' @param size Point size.
#' @param width Width of the figure.
#' @param height Height of the figure.
#' @param rmelements A vector of elememts to be removed {\link[ggpubr]{rremove}}.
#' @param ... Additional arguments passed on to \code{\link[SOT]{umap}}
#' @return 
#' \itemize{
#'   \item Coordinate of 2-D layout
#'   \item ggplot2 object
#'   \item igraph 
#' }
#' @export
layout_2d <- function(sce,
                   usedEmbed = "PCA",
                   usedDims = NULL,
                   datatype = NULL,
                   filename = NULL,
                   k = 35,
                   n_trees = 50,
                   max_iter = 2,
                   seed = 1, 
                   layout_method = "FR",
                   grid = "auto",
                   iterations = 1000,
                   repel = 2,
                   plotstep = 0,
                   clus_method = NULL,
                   shape = NULL,
                   colour_condition = NULL,
                   colour_gene = NULL,
                   colour_embed = NULL,
                   color_scale = NULL,
                   cmap = NULL, 
                   compute_layout = TRUE,
                   orderExp = FALSE,
                   title = NULL,
                   ncol = NULL, 
                   nrow = NULL,
                   alpha = 0.9,
                   size = 1,
                   width = 4,
                   height = 3,
                   rmelements = NULL,
                   ...){
  if (class(sce) != "SingleCellExperiment"){
    stop("sce must be SingleCellExperiment object")
  }
  else if (!usedEmbed %in% reducedDimNames(sce)){
    stop(usedEmbed, "is not in reduceDims slot, current available embedding are", reducedDimNames(sce))
  }
  else{
    mat <- reducedDim(sce, usedEmbed)
  }
  if (!is.null(usedDims)){
    if (!all(unique(usedDims) %in% 1:ncol(mat))){
      stop("usedDims are not matched the column number of usedEmbed")
    }
    else{
      mat <- mat[, usedDims]
    }
  }
  else{
    usedDims <- 1:ncol(mat)
  }
  if (is.null(seed)){
    seed <- 1
  }
  suffix = ""
  if (!is.null(colour_gene)){
    if (is.numeric(colour_gene)){
      if (length(colour_gene) == nrow(mat)){
        expression = data.frame(Level = colour_gene)
        genelv <- "Value"
      }
    }
    else if (!any(colour_gene %in% rownames(sce))){
      stop("Gene specified by colour_gene is not in sce object")
    }
    else{
      colour_gene <- colour_gene[colour_gene %in% rownames(sce)]
      genelv <- unique(colour_gene)
      # assay(sce, datatype) <- as(assay(sce, datatype), "matrix")
      expression <- as.data.frame(t(as.matrix(assay(sce[colour_gene, ], datatype))))
    }
    suffix <- ".gene"
  }
  
  if (!is.null(colour_embed)){
    stopifnot(any(colour_embed %in% reducedDimNames(sce)))
    stopifnot(length(colour_embed) == 1)
    expression <- as.data.frame(reducedDim(sce, colour_embed))
    genelv <- sort(colnames(expression))
    cmap <- "Spectral"
    suffix <- paste0(".",colour_embed)
  }
  

  # Construct knn graph
  knn_idx <- knn.index(mat, k = k)
  al <- split(knn_idx, 1:nrow(knn_idx))
  g <- graph_from_adj_list(al, mode = "out")
  colData(sce)$degree <- degree(g)
  
  # Clustering
  # TODO: Add more clustering methods
  if (!is.null(clus_method)){
    if (clus_method == "louvain"){
      message(clus_method, " clustering on ", "kNN graph (k = ", k, ")...")
      gcom <- cluster_louvain(as.undirected(g))
    }
    else if (clus_method == "spinglass"){
      message(clus_method, " clustering on ", "kNN graph (k = ", k, ")...")
      gcom <- cluster_spinglass(as.undirected(g))
    }
    else{
      stop("Availabel clustering methods are ", paste("<louvain,", "spinglass>", sep = " "))
    }
    if (is.null(colour_condition)){
      colour_condition <- clus_method
    }
    colData(sce)[, clus_method] <- factor(gcom$membership)
    message("Finished! ", length(unique(gcom$membership)), " communities found.")
  }
  
  if (!is.null(colour_condition)){
    if (!colour_condition %in% colnames(colData(sce))){
      stop("colour_condition is not in colData. The avalable conditions are", colnames(colData(sce)))
    }
    else{
      condition <- colData(sce)[, colour_condition]
    }
    if (!is.null(color_scale)){
      if (length(unique(color_scale)) != length(unique(condition))){
        stop("color_scale and are not matched")
      }
    }
    else if (!is.null(cmap)){
      color_scale <- cmap
    }
    else{
      #color_scale <- brewer.pal(n = 12, name = "Paired")
      color_scale <- brewer.pal(n = 8, name = "Accent")
    }
  }
  
  # Graph layout
  if (any(grep("FR|largeVis|umap", usedEmbed)) | !compute_layout){
    xy = as.data.frame(reducedDim(sce, usedEmbed))
    rd.name <- usedEmbed
  }
  else{
    set.seed(seed)
    if (layout_method == "FR"){
      xy <- layout_with_fr(g, grid = grid, dim = 2, niter = iterations)
    }
    else if (layout_method == "largeVis"){
      if (!"largeVis" %in% rownames(installed.packages())){
        if (!"devtools" %in% rownames(installed.packages())){
          install.packages("devtools")
        }
        devtools::install_github("analyxcompany/ForceAtlas2")
      }
      xy <- t(largeVis::largeVis(t(reducedDim(sce, usedEmbed)), K=k, n_trees=n_trees, max_iter=max_iter, seed=seed)$coords)
    }
    else if (layout_method == "umap"){
      xy <- umap(reducedDim(sce, usedEmbed), random_state=seed, ...)
    }
    else{
      stop("Layout method should be one of FR, largeVis or umap")
    }
    xy <- as.data.frame(xy)
    rd.name <- paste0(usedEmbed, "_", layout_method)
  }
  colnames(xy) <- c("Dim 1", "Dim 2")

  # Set plot style
  mytheme <- theme_bw() + theme(legend.position="right",
                                panel.grid.minor = element_blank(),
                                panel.grid.major = element_blank(),
                                axis.text = element_blank(),
                                axis.ticks = element_blank()) # panel.border = element_blank() 
  
  if (!is.null(colour_gene) | !is.null(colour_embed)){
    if (is.null(cmap)){
      colour_gradientn <- c("#999999","#999999","#e60000","#e60000")
    }
    else{
      colour_gradientn <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(255)
    }
    # GeneExp <- cbind(xy, expression) %>%
    #   gather(Gene, Expression, -`Dim 1`, -`Dim 2`)
    # GeneExp$Gene <- factor(GeneExp$Gene, levels = genelv)
    plist = list()
    for (vs in colnames(expression)){
      xy$value = expression[, vs]
      if (orderExp){
        xy <- xy[order(xy$value, decreasing = F), ]
      }
      p0 <- ggplot(xy, aes(`Dim 1`, `Dim 2`)) + 
        geom_point(aes(colour = value), alpha = alpha, size = size, stroke = 0) + 
        labs(title = vs, color = NULL, x = "Dim 1", y = "Dim 2") + 
        scale_colour_gradientn(colours = colour_gradientn) +
        coord_fixed(ratio = (max(xy$`Dim 1`)-min(xy$`Dim 1`))/(max(xy$`Dim 2`)-min(xy$`Dim 2`))) + mytheme
      if (!is.null(rmelements)){
        for (el in rmelements){
          p0 <- p0 + rremove(el) 
        }
      }
      plist[[vs]] <- p0
    }
    p <- ggarrange(plotlist = plist, ncol = ncol, nrow = nrow)
  }
  else if (!is.null(colour_condition)){
    suffix <- paste0(".",colour_condition)
    p <- ggplot(xy, aes(`Dim 1`, `Dim 2`)) + 
         geom_point(aes(colour = condition), alpha = alpha, size = size, stroke = 0) + 
         labs(color = colour_condition,x = "Dim 1", y = "Dim 2") +
         scale_colour_manual(values = colorRampPalette(color_scale)(length(unique(condition)))) +
         coord_fixed(ratio = (max(xy$`Dim 1`)-min(xy$`Dim 1`))/(max(xy$`Dim 2`)-min(xy$`Dim 2`))) + 
        guides(colour = guide_legend(override.aes = list(size=3))) + mytheme
  }
  else{
    p <- ggplot(xy, aes(`Dim 1`, `Dim 2`)) + 
         geom_point(alpha = alpha, size = size, stroke = 0) + 
         labs(color = "None", x = "Dim 1", y = "Dim 2") +
         coord_fixed(ratio = (max(xy$`Dim 1`)-min(xy$`Dim 1`))/(max(xy$`Dim 2`)-min(xy$`Dim 2`))) +
         guides(colour = guide_legend(override.aes = list(size=3))) + mytheme
  }
  reducedDim(sce, rd.name) <- as.matrix(xy[, c("Dim 1", "Dim 2")])
  metadata(sce)[[paste0("p.",rd.name,suffix)]] <- p
  message("Use metadata(sce)$",paste0("p.",rd.name,suffix)," to visualize embeddings.")
  if (!is.null(filename)){
    ggsave(filename, width = width, height = height)
  }
  return(sce)
}

#' Heatmap of the reduced dimension data
#' 
#' Heatmap of the reduced dimension data
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom pheatmap pheatmap
#' @importFrom fastcluster hclust
#' @param sce A SingleCellExperiment object.
#' @param usedEmbed Specify the reduce dimension type for construct network. See current avaible embedding by \code{\link[SingleCellExperiment]{reducedDimNames}}.
#' @param k Number of clusters.
#' @param method Clustering method.
#' @param filename Figure name.
#' @param anno_col Characters to specify annotate conditions.
#' @param color_scale Customer color scale for each condition.
#' @param cmap Color map.
#' @return Clustering result.
#' @export
plot_embed <- function(sce,
                       usedEmbed,
                       k = 8,
                       method = "complete",
                       filename = NA,
                       cells.use = NULL,
                       vars.use = NULL,
                       anno_col = NULL,
                       color_scale = NULL,
                       cmap = NULL,
                       cluster_row = TRUE,
                       ...){
  rdm <- as.matrix(reducedDim(sce, usedEmbed))
  if (!is.null(cells.use)){
    rdm <- rdm[cells.use,]
  }
  if (!is.null(vars.use)){
    rdm <- rdm[,vars.use]
  }
  # Clustering
  hc <- hclust(dist(rdm), method = method)
  memb <- as.data.frame(factor(cutree(hc, k = k)))
  colnames(memb) <- "Cluster"
  co <- hc$order
  if (is.null(color_scale)){
    if (!is.null(cmap)){
      n = gsub("*.([0-9]*)", "\\1", cmap)
      name = strsplit(cmap, n)[[1]]
      color = brewer.pal(n = as.numeric(n), name = name)
    }
    else{
      color = brewer.pal(n = 12, name = "Paired")
    }
  }
  else{
    color <- color_scale
  }
  annotation_colors <- list()
  if (is.null(anno_col)){
    condition <- memb
    annotation_colors[[colnames(memb)]] <- colorRampPalette(color)(length(unique(memb[, 1])))
    names(annotation_colors[[colnames(memb)]]) <- levels(factor(memb[, 1]))
  }
  else{
    stopifnot(any(anno_col %in% colnames(colData(sce))))
    anno_col <- anno_col[anno_col %in% colnames(colData(sce))]
    condition <- colData(sce)[, anno_col]
    condition <- cbind(condition, memb)
    for (cond in colnames(condition)){
      annotation_colors[[cond]] <- colorRampPalette(color)(length(unique(condition[, cond])))
      names(annotation_colors[[cond]]) <- levels(factor(condition[, cond]))
    }
  }

  pheatmap(t(rdm[co, ]), 
           filename = filename,
           cluster_rows = cluster_row,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           annotation_col = condition,
           annotation_colors = annotation_colors,
           ...)
  
  return(invisible(memb))
  
}

dlm_helper <- function(g){
  t <- as.ts(as.numeric(g), start = 1, end = length(g))
  nileBuild <- function(par) {
    dlmModPoly(1, dV = exp(par[1]), dW = exp(par[2]))
  }
  nileMLE <- dlmMLE(t, rep(0,2), nileBuild)
  nileMod <- nileBuild(nileMLE$par)
  nileFilt <- dlmFilter(t, nileMod)
  nileSmooth <- dlmSmooth(nileFilt)
  
  return(nileSmooth$s[-1])
  
}

#' Uniform Manifold Approximation and Projection (UMAP).
#' 
#' Dimension reduction technique that can be used for visualisation similarly to t-SNE, but also for general non-linear dimension reduction.
#' @importFrom reticulate py_module_available import
#' @param mat Data with rows are samples and cols are genes.
#' @param n_neighbors float (optional, default 15) The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100.
#' @param n_components int (optional, default 2) The dimension of the space to embed into.
#' @param metric string (optional, default 'euclidean') The metric to use to compute distances in high dimensional space. Valid string metrics include: euclidean, manhattan, minkowski, mahalanobis, seuclidean, cosine, correlation, hamming, jaccard, kulsinski.
#' @param negative_sample_rate int (optional, default 5) The number of negative edge/1-simplex samples to use per positive edge/1-simplex sample in optimizing the low dimensional embedding.
#' @param init string (optional, default 'spectral') How to initialize the low dimensional embedding. Options are: 'spectral', 'random', A numpy array of initial embedding positions.
#' @param min_dist float (optional, default 0.1) The effective minimum distance between embedded points. The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out.
#' @param random_state int, RandomState instance or NULL, optional (default: NULL)
#' @param verbose bool (optional, default False) Controls verbosity of logging.
#' @return Low dimension embedding
#' @export
umap = function(mat, 
                umap.args = list(),
                n_neighbors = 20,
                n_components = 2, 
                metric = "euclidean",
                negative_sample_rate = 5,
                init = "spectral",
                min_dist = 0.3,
                random_state = 1,
                verbose = FALSE){
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through conda or pip (e.g. conda install -c conda-forge umap-learn).")
  }
  message("UMAP for ", nrow(mat), " samples and ", ncol(mat), " features, embedding in ", n_components, " dimensions.")
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = n_neighbors),
    n_components = as.integer(x = n_components),
    metric = metric,
    negative_sample_rate = negative_sample_rate,
    init = init,
    min_dist = min_dist,
    random_state = as.integer(x = random_state)
  )
  embedding <- umap$fit_transform(as.matrix(x = mat))
  colnames(embedding) <- paste0("UMAP", 1:n_components)
  rownames(embedding) <- rownames(mat)
  
  return(embedding)
  
}
