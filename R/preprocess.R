#' Filtering low expression
#' 
#' Filtering the genes expressed at least minexp in less than mincells.
#' @importFrom SummarizedExperiment assay assays rowData
#' @import dplyr
#' @param sce SingleCellExperiment object or a matrix.
#' @param datatype Sepcify the data type in sce for filtering.
#' @param minexp Minimum expression level.
#' @param mincells Minimum cells expressed minexp.
#' @return A sce object of a matrix aftering filtering.
#' @export
FilterLow = function(sce, datatype = NULL, minexp = 10, mincells = 5){
  if (!class(sce) == "SingleCellExperiment"){
    stop("sce should be a SingleCellExperiment object")
  }
  if (!is.null(datatype)){
    if (!datatype %in% names(assays(sce))){
      stop("Available datatype are", names(assays(sce)))
    }
  }
  data = assay(sce, datatype)
  fillow.mask <- rowSums(data >= minexp) >= mincells
  rowData(sce)[, "filter.low"] <- fillow.mask
  if ("genes.use" %in% colnames(rowData(sce))){
    rowData(sce)[, "genes.use"] <- rowData(sce)[, "genes.use"] & fillow.mask
  }else{
    rowData(sce)[, "genes.use"] <- fillow.mask
  }
  if (!"symbol" %in% colnames(rowData(sce))){
    rowData(sce) <- cbind(data.frame(symbol=rownames(sce)), as.data.frame(rowData(sce)))
  }
  rowData(sce) <- as.data.frame(rowData(sce)) %>% select(-genes.use, genes.use)
  return(sce)
}

#' Aderson-Darling test
#' 
#' Test if gene expression in different sampling time points come from a common population
#' @importFrom SummarizedExperiment assay assays rowData
#' @importFrom kSamples ad.test
#' @importFrom parallel makeCluster
#' @importFrom stats p.adjust
#' @importFrom reticulate source_python
#' @importFrom S4Vectors metadata
#' @import doParallel
#' @import foreach
#' @import dplyr
#' @param sce SingleCellExperiment object.
#' @param datatype Sepcify the data type in sce for filtering.
#' @param condition Conditions corresponding to cells.
#' @param useLevels Subset conditions for test - default is NULL.
#' @param genes.use A logical vector in rowData to specify the genes for calculation.
#' @param sample.cells Sampling cell number in each condition to speed up testing.
#' @param adj.method p-value correction method. See \code{\link[stats]{p.adjust}} for more detail.
#' @param r.implement Whether to use R implement of AD test. If FALSE, I will use anderson_ksamp function in scipy.
#' @param ncore Number of cores used for parallel. 
#' @return Adjusted p-value of Anderson-Darling test.
#' @export
ADtest <- function(sce, 
                   datatype = "logcounts",
                   condition, 
                   useLevels = NULL,
                   genes.use = "find.hvgs",
                   sample.cells = 500,
                   adj.method = "BH",
                   thr.padj = 1e-3,
                   r.implement = FALSE,
                   ncore = 6){
  if (class(sce) != "SingleCellExperiment"){
      stop("sce must be a sce object")
  }
  if (!is.null(datatype)){
    if (!datatype %in% names(assays(sce))){
      stop("Available datatype are", names(assays(sce)))
    }
  }
  if (!condition %in% colnames(colData(sce))){
    stop("cond must be in colData")
  }
  else{
    cond = as.vector(colData(sce)[, condition])
  }
  if (!is.null(useLevels)){
    if (!all(useLevels %in% cond)){
      stop("useLevels must be in ", paste(unique(cond), sep = " ", collapse = T))
    }
  }
  else{
    useLevels <- unique(as.vector(cond))
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
  if (!"symbol" %in% colnames(rowData(sce))){
    rowData(sce) <- cbind(data.frame(symbol=rownames(sce)), as.data.frame(rowData(sce)))
  }
  lv <- cond[cond %in% useLevels]

  # Sampling cells to speed up testing
  cells.use <- colnames(sce)[cond %in% useLevels]
  names(lv) <- cells.use

  if (!is.null(sample.cells)){
    set.seed(1)
    cells.use <- do.call(c, lapply(useLevels, function(i) {sample(cells.use[lv == i], min(sum(lv == i), sample.cells))}))
    lv <- lv[cells.use]
  }
  
  vi <- apply(assay(sce[, cells.use], datatype), 1, function(x) sum(abs(diff(x))) >0 )
  genes.use <- genes.use & vi
  message("Perform Aderson-Darling test for ", length(cells.use), " cells of <", paste(useLevels, collapse = " "),">")
  if (ncore == 1){
    message("Use 1 core to perform Anderson-Darling test...")
    if (!py_module_available(module = "scipy")) {
      pv = apply(assay(sce[genes.use, cells.use], datatype), 1, function(y) ad.test(y ~ as.character(lv))$ad[1,3]) # p-value of asymptotic method
    }else{
      # The scipy application is much faster
      message("*Use python implement of anderson_ksamp*")
      scipy <- import(module = "scipy", delay_load = TRUE)
      pv <- apply(assay(sce[genes.use, cells.use], datatype), 1, function(y) {
        l <- split(as.numeric(y), as.character(lv))
        names(l) <- NULL
        scipy$stats$anderson_ksamp(l)$significance_level
      })
      }
  }
  else{
    message("Use ", ncore ," cores to perform Anderson-Darling test...")
    cl = makeCluster(ncore)
    registerDoParallel(cl)
    if (r.implement) {
      pv <- unlist(foreach(i = 1:sum(genes.use), .packages=c("kSamples", "SummarizedExperiment")) %dopar% {
        ad.test(as.numeric(assay(sce[genes.use, cells.use][i, ], datatype)) ~ as.character(lv))$ad[1,3] # p-value of asymptotic method
      })
      
    }else{
      if (py_module_available(module = "scipy")){
        message("*Use python implement of anderson_ksamp*")
        idx <- 1:sum(genes.use)
        block <- as.numeric(cut(idx, breaks = ncore))
        pv <- unlist(foreach(i = 1:ncore, .packages=c("reticulate", "SummarizedExperiment"), .combine=c) %dopar% {
          source_python(system.file("py_script", "adtest.py", package = "SOT"))
          py_adtest(assay(sce[genes.use, cells.use][block == i, ], datatype), as.character(lv))
          # l <- split(as.numeric(assay(sce[genes.use, cells.use][i, ], datatype)), as.character(lv))
          # names(l) <- NULL
          # scipy$stats$anderson_ksamp(l)$significance_level
        })
      }else{
        stop("Python package scipy is not available.")
      }
    }
    stopCluster(cl)
  }
  adtest.padj <- rep(NA, nrow(sce))
  names(adtest.padj) <- rownames(sce)
  adtest.padj[rownames(sce[genes.use, ])] <- p.adjust(pv, method = adj.method)
    
  rowData(sce)$`adtest.padj` <- adtest.padj
  metadata(sce)$`AD test condistions` <- useLevels
  
  sig.genes <- adtest.padj[adtest.padj < thr.padj]
  sig.genes <- sig.genes[!is.na(sig.genes)]
  ad.mask <- rownames(sce) %in% names(sig.genes)
  
  rowData(sce)[, "ad.test"] <- ad.mask
  if ("genes.use" %in% colnames(rowData(sce))){
    rowData(sce)[, "genes.use"] <- rowData(sce)[, "genes.use"] & ad.mask
  }else{
    rowData(sce)[, "genes.use"] <- ad.mask
  }
  rowData(sce) <- as.data.frame(rowData(sce)) %>% select(-genes.use, genes.use)
  return(sce)
}

#' Find high variable genes
#' 
#' Fing high variable genes by fitting a mean-dependent trend to the gene-specific variances.
#' @references \url{https://f1000research.com/articles/5-2122/v2}
#' @importFrom scran trendVar decomposeVar
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom graphics curve points
#' @import dplyr
#' @param sce SingleCellExperiment object.
#' @param datatype Sepcify the data type in sce to find hvgs.
#' @param genes.use A logical vector in rowData to specify the genes for calculation.
#' @param thr.bio The threshold of biological component of the variance.
#' @param thr.tech The threshold of Technical component of the variance.
#' @param thr.p.value The threshold of p-values for the test against the null hypothesis that bio=0.
#' @param thr.FDR The threshold of adjusted p-values.
#' @param thr.low The low threshold of mean log-expression.
#' @param thr.high The high threshold of mean log-expression.
#' @param show.plot Whether to plot the fitting results.
#' @param ... Additional arguments passed on to \code{\link[scran]{trendVar}}.
#' @return SingleCellExperiment object with hvgs mask in rowData slot.
#' @export
FindHVGs <- function(sce, 
                     datatype = "logcounts", 
                     genes.use = "filter.low", 
                     thr.bio = 0,
                     thr.tech = NULL,
                     thr.p.value = NULL,
                     thr.FDR = 0.1,
                     thr.low = 0,
                     thr.high = 20,
                     show.plot = TRUE,
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
  if (!"symbol" %in% colnames(rowData(sce))){
    rowData(sce) <- cbind(data.frame(symbol=rownames(sce)), as.data.frame(rowData(sce)))
  }
  var.fit <- trendVar(assay(sce[genes.use, ], datatype), ...) # loess.args=list(span=0.3), parametric=T
  var.out <- decomposeVar(assay(sce[genes.use, ], datatype), var.fit)
  hvg <- var.out
  if (!is.null(thr.bio)){
    hvg = hvg[hvg$bio > thr.bio,]
  }
  if (!is.null(thr.tech)){
    hvg <- hvg[hvg$tech < thr.tech,]
  }
  if (!is.null(thr.p.value)){
    hvg <- hvg[hvg$p.value < thr.p.value,]
  }
  if (!is.null(thr.FDR)){
    hvg <- hvg[hvg$FDR < thr.FDR,]
  }
  if (!is.null(thr.low)){
    hvg <- hvg[hvg$mean > thr.low,]
  }
  if (!is.null(thr.high)){
    hvg <- hvg[hvg$mean < thr.high,]
  }
  if (show.plot){
    plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
         ylab="Variance of log-expression", main = "Fit a variance trend")
    curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
    points(hvg$mean, hvg$total, col = "red", pch = 16, cex=0.6)
  }
  hvg.mask <- rownames(sce) %in% rownames(hvg)
  rowData(sce)[, "find.hvgs"] <- hvg.mask
  if ("genes.use" %in% colnames(rowData(sce))){
    rowData(sce)[, "genes.use"] <- rowData(sce)[, "genes.use"] & hvg.mask
  }else{
    rowData(sce)[, "genes.use"] <- hvg.mask
  }
  fit.df <- data.frame(matrix(rep(NA, ncol(var.out)*nrow(sce)),ncol = ncol(var.out)),row.names = rownames(sce))
  colnames(fit.df) <- colnames(var.out)
  fit.df[rownames(var.out), ] <- as.data.frame(var.out)
  # rowData(sce) <- cbind(rowData(sce), fit.df[rownames(sce), ])
  rowData(sce) <- as.data.frame(rowData(sce)) %>% select(-genes.use, genes.use)
  return(sce)
}
