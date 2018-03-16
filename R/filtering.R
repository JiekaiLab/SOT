#' Filtering low expression
#' 
#' Filtering the genes expressed at least minexp in less than mincells.
#' @importFrom SummarizedExperiment assay assays
#' @param sce SingleCellExperiment object or a matrix.
#' @param minexp Minimum expression level.
#' @param mincells Minimum cells expressed minexp.
#' @param datatype Sepcify the data type in sce for filtering.
#' @return A sce object of a matrix aftering filtering.
#' @export
FilterLow = function(sce, minexp = 10, mincells = 5, datatype = NULL){
  if (class(sce) == "SingleCellExperiment"){
    if (!is.null(datatype)){
      if (!datatype %in% names(assays(sce))){
        stop("Available datatype are", names(assays(sce)))
      }
    }
    data = assay(sce, datatype)
    sce = sce[rowSums(data >= minexp) >= mincells, ]
    
    return(sce)
    
  }
  else{
    return(sce[rowSums(sce >= minexp) >= mincells, ])
  }
}

#' Aderson-Darling test
#' 
#' Test if gene expression in different sampling time points come from a common population
#' @importFrom kSamples ad.test
#' @importFrom parallel makeCluster
#' @importFrom stats p.adjust
#' @import doParallel
#' @import foreach
#' @param sce SingleCellExperiment object.
#' @param condition Conditions corresponding to cells.
#' @param useLevels Subset conditions for test - default is NULL.
#' @param adj.method p-value correction method. See \code{\link[stats]{p.adjust}} for more detail.
#' @param datatype Sepcify the data type in sce for filtering.
#' @param ncore Number of cores used for parallel. 
#' @return Adjusted p-value of Anderson-Darling test.
#' @export
adtest <- function(sce, condition, useLevels = NULL, adj.method = "BH", datatype = NULL, ncore = 6){
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
    cond = colData(sce)[, condition] 
  }
  if (!is.null(useLevels)){
    if (!all(useLevels %in% cond)){
      stop("useLevels must be in ", paste(unique(cond), sep = " ", collapse = T))
    }
  }
  else{
    useLevels <- unique(cond)
  }
  lv <- cond[cond %in% useLevels]
  vi <- which(abs(apply(assay(sce[, cond %in% useLevels], datatype), 1, function(x) sum(abs(diff(x))))) > 0)
  sub_sce <- sce[vi, ]
  
  message("Perform Aderson-Darling test for <", paste(useLevels, collapse = " "),">")
  if (ncore == 1){
    pv = apply(assay(sub_sce[, cond %in% useLevels], datatype), 1, function(y) ad.test(y ~ as.character(lv))$ad[1,3]) # p-value of asymptotic method
  }
  else{
    cl = makeCluster(ncore)
    registerDoParallel(cl)
    pv <- unlist(foreach(i = 1:nrow(assay(sub_sce, datatype))) %dopar% {
      ad.test(as.numeric(assay(sub_sce[i, cond %in% useLevels], datatype)) ~ as.character(lv))$ad[1,3] # p-value of asymptotic method
    })
  }
  rowData(sub_sce)$`adtest.padj` = p.adjust(pv, method = adj.method)
  metadata(sub_sce)$`AD test condistions` = useLevels
  
  return(sub_sce)
  
}
