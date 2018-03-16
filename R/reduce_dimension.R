#' Force-directed layout
#' 
#' 2-D visualization of data on kNN graph using force-directed layout
#' @importFrom FNN knn.index
#' @importFrom igraph graph_from_adj_list layout_with_fr 
#' @importFrom stats cor prcomp
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SingleCellExperiment reducedDim reducedDimNames colData
#' @importFrom SummarizedExperiment colData
#' @import ggplot2
#' @import tidyr
#' @param sce A SingleCellExperiment object.
#' @param usedEmbed Specify the reduce dimension type for construct network.
#' @param usedDims A vector to specify particular dimension in usedEmbed to calculate.
#' @param datatype Type of expression.
#' @param filename Filename to save - default is NULL.
#' @param k k nearest neighbors to construct graph.
#' @param seed Random seed.
#' @param colour_condition Label the cells with continuous color.
#' @param colour_gene A character of a vector of a gene expression.
#' @param color_scale Scatter plot color correponding to levels - default is NULL.
#' @param cmap Color style of cmap will be used when color_scale is NULL.
#' @param title Plot title.
#' @param alpha Transparency of points.
#' @param size Point size.
#' @param width Width of the figure.
#' @param height Height of the figure.
#' @return 
#' \itemize{
#'   \item Coordinate of 2-D layout
#'   \item ggplot2 object
#' }
#' @export
fdl_2d <- function(sce,
                   usedEmbed = "PCA",
                   usedDims = NULL,
                   datatype = NULL,
                   filename = NULL,
                   k = 30,
                   seed = 1, 
                   colour_condition = NULL,
                   colour_gene = NULL,
                   color_scale = NULL,
                   cmap = "Paired", 
                   title = "Force-directed layout",
                   alpha = 0.9,
                   size = 1,
                   width = 4,
                   height = 3){
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
  if (!is.null(colour_gene)){
    if (is.numeric(colour_gene)){
      if (length(colour_gene) == nrow(mat)){
        expression = data.frame(Level = colour_gene)
        genelv <- "Level"
      }
    }
    else if (!any(colour_gene %in% rownames(sce))){
      stop("Gene specified by colour_gene is not in sce object")
    }
    else{
      colour_gene <- colour_gene[colour_gene %in% rownames(sce)]
      genelv <- unique(colour_gene)
      expression <- as.data.frame(t(assay(sce[colour_gene, ], datatype)))
    }
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
      color_scale <- brewer.pal(n = length(unique(condition)), name = cmap)
    }
    else{
      color_scale <- NULL # Is there better choice for color?
    }
  }

  # Construct knn graph
  knn_idx <- knn.index(mat, k = k)
  al <- split(knn_idx, 1:nrow(knn_idx))
  g <- graph_from_adj_list(al, mode="out")
  
  # TODO: Add more layout methods
  # Graph layout
  if (any(grep("FR", usedEmbed))){
    xy = as.data.frame(reducedDim(sce, usedEmbed))
  }
  else{
    set.seed(seed)
    co <- layout_with_fr(g, grid = "nogrid", dim = 2)
    xy <- as.data.frame(co)
    colnames(xy) <- c("Dim 1", "Dim 2")
  }

  # Set plot style
  mytheme <- theme_bw() + theme(legend.position="bottom",
                                panel.grid.minor = element_blank(),
                                panel.grid.major = element_blank(),
                                axis.text = element_blank(),
                                axis.ticks = element_blank()) 
  
  if (!is.null(colour_gene)){
    GeneExp <- cbind(xy, expression) %>%
      gather(Gene, Expression, -`Dim 1`, -`Dim 2`)
    GeneExp$Gene <- factor(GeneExp$Gene, levels = genelv)
    p <- ggplot(GeneExp, aes(`Dim 1`, `Dim 2`)) + 
         facet_wrap(~ Gene) +
         geom_point(aes(colour = Expression), alpha = alpha, size = size) + 
         labs(color = "Expression level", title = title) + 
         scale_colour_gradientn(colours = c("#999999","#999999","#e60000","#e60000")) +
         coord_fixed(ratio = (max(xy$`Dim 1`)-min(xy$`Dim 1`))/(max(xy$`Dim 2`)-min(xy$`Dim 2`)))
         # scale_colour_gradientn(colors = myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(255)) 
  }
  else if (!is.null(colour_condition)){
    p <- ggplot(xy, aes(`Dim 1`, `Dim 2`)) + 
         geom_point(aes(colour = condition), alpha = alpha, size = size) + 
         labs(color = "Day", title = title) +
         scale_colour_manual(values = color_scale) +
         coord_fixed(ratio = (max(xy$`Dim 1`)-min(xy$`Dim 1`))/(max(xy$`Dim 2`)-min(xy$`Dim 2`)))
  }
  else{
    p <- ggplot(xy, aes(`Dim 1`, `Dim 2`)) + 
         geom_point(alpha = alpha, size = size) + 
         labs(color = "None", title = title) +
         coord_fixed(ratio = (max(xy$`Dim 1`)-min(xy$`Dim 1`))/(max(xy$`Dim 2`)-min(xy$`Dim 2`)))
  }
  p <- p + mytheme
  reducedDim(sce, paste("FR", usedEmbed, paste(usedDims, collapse = ""), sep = "_")) <- as.matrix(xy)
  if (!is.null(filename)){
    ggsave(filename, width = width, height = height)
  }
  
  return(list(sce = sce, p = p))
  
}
