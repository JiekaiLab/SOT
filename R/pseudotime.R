#' Wishbone pseudotime
#' 
#' Graph based pasudotime infer. Please install python package "wishbone" before using this function.
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment reducedDim
#' @param sce A SingleCellExperiment object or matrix.
#' @param usedEmbed Specify the reduce dimension type for construct network. See current avaible embedding by \code{\link[SingleCellExperiment]{reducedDimNames}}.
#' @param start_cell Character specified a start cell of the trajectory.
#' @param datatype Type of expression used for plot. See current avaible assays by \code{\link[SummarizedExperiment]{assays}}.
#' @param usedDims A vector to specify particular dimensions in usedEmbed to calculate.
#' @param condition Used with usedLevels. Specify a condition in colData(sce) for trimming cells.
#' @param usedLevels Used with condition. Only the cells in usedLevels are used for calculating.
#' @param plot_genes A vector of genes id. If plot_genes is not NULL, plot the average smoothing gene expression along the trajectories.
#' @param branch Whether branches exist or not.
#' @param k kNN graph.
#' @param num_waypoints Number of num_waypoints.
#' @param seed Random seed.
#' @param figsize A vector correponding to the figure width and height.
#' @return sce object contain pseudotime and branches.
#' @export
wishbone <- function(sce,
                     usedEmbed,
                     start_cell, 
                     datatype = "scaledata",
                     usedDims = NULL,
                     condition = NULL, 
                     usedLevels = NULL,
                     plot_genes = NULL,
                     branch = TRUE,
                     k = 35,
                     num_waypoints = 150,
                     seed = 1,
                     figsize = c(20, 10)){
  if (!is.null(usedLevels)){
    if (!is.null(condition)){
      sce_trim <- sce[, colData(sce)[, condition] %in% usedLevels]
    }
  }
  else{
    sce_trim <- sce
  }
  if (!start_cell %in% colnames(sce_trim)){
    stop("start_cell is not in the list")
  }

  exprmat <- t(assay(sce_trim, datatype))
  rdmat <- reducedDim(sce_trim, usedEmbed)
  if (!is.null(usedDims)){
    if (max(usedDims) <= ncol(rdmat)){
      rPython::python.assign("components_list", usedDims-1)
    }
    else{
      stop("The maximum dimensions of usedDims is ", ncol(rdmat))
    }
  }
  else{
    rPython::python.assign("components_list", 1:ncol(rdmat)-1)
  }
  if (branch){
    rPython::python.exec("branch=bool(1)")
  }
  else{
    rPython::python.exec("branch=bool(0)")
  }
  message("Running wishbone ...")
  rPython::python.exec("import numpy as np;import wishbone;import matplotlib.pyplot as plt;import pandas as pd")
  rPython::python.assign("genes", colnames(exprmat))
  rPython::python.assign("cells", rownames(exprmat))
  rPython::python.assign("dimnames", colnames(rdmat))
  rPython::python.assign("exprmat", as.vector(exprmat))
  rPython::python.assign("rdmat", as.vector(rdmat))
  rPython::python.assign("start_cell", start_cell)
  rPython::python.assign("rcdim", dim(exprmat))
  rPython::python.assign("rddim", dim(rdmat))
  rPython::python.assign("k", k)
  rPython::python.assign("num_waypoints", num_waypoints)
  rPython::python.assign("seed", seed)
  rPython::python.exec("exprmat=np.array(exprmat);exprmat=np.reshape(exprmat,rcdim,order='F');exprmat=pd.DataFrame(data=exprmat, index=cells, columns=genes)")
  rPython::python.exec("rdmat=np.array(rdmat);rdmat=np.reshape(rdmat,rddim,order='F');rdmat=pd.DataFrame(data=rdmat, index=cells, columns=dimnames)")
  rPython::python.exec("scdata=wishbone.wb.SCData(exprmat);scdata.diffusion_eigenvectors=rdmat;np.random.seed(seed);wb=wishbone.wb.Wishbone(scdata);wb.run_wishbone(start_cell=start_cell, branch=branch, components_list=components_list, k=k, num_waypoints=num_waypoints)")
  tra <- rPython::python.get("list(wb.trajectory)")
  br <- rPython::python.get("list(wb.branch)")

  if (!is.null(plot_genes)){
    plot_genes <- plot_genes[plot_genes %in% rownames(sce_trim)]
    if (length(plot_genes) > 0){
      rPython::python.assign("plot_genes", plot_genes)
      rPython::python.assign("w", figsize[1])
      rPython::python.assign("h", figsize[2])
      rPython::python.exec("vals, fig, ax = wb.plot_marker_trajectory(plot_genes);fig.set_size_inches(w, h);plt.savefig('wishbone_plot_genes.pdf')")
    }
  }
  # Store results
  wb_df <- data.frame(trajectory = rep(0, ncol(sce)), branch = rep(0, ncol(sce)), row.names = colnames(sce))
  wb_df[colnames(sce_trim), ] <- cbind(tra, br)
  colData(sce)$trajectory <- wb_df$trajectory
  colData(sce)$branch <- factor(wb_df$branch)
  
  return(sce)
  
}

#' Bayesian and Likelihood Analysis of Dynamic Linear Models
#' 
#' @importFrom dlm dlmModPoly dlmMLE dlmFilter dlmSmooth
dlm_helper <- function(g){
  t <- as.ts(as.numeric(g), start = 1, end = length(g))
  exprBuild <- function(par) {
    dlmModPoly(1, dV = exp(par[1]), dW = exp(par[2]))
  }
  exprMLE <- dlmMLE(t, rep(0,2), exprBuild)
  exprMod <- exprBuild(exprMLE$par)
  exprFilt <- dlmFilter(t, exprMod)
  exprSmooth <- dlmSmooth(exprFilt)
  
  return(exprSmooth$s[-1])
  
}

#' Smoothing expression
#' 
#' Smoothing expression using Bayesian and Likelihood Analysis of Dynamic Linear Models.
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @param mat Expression matrix.
#' @param bincells Use the mean of bincells for smoothing.
#' @param normalize Normalize the smoothed value to 0-1.
#' @param ncore ncores for parallel.
#' @return Smoothing expression.
#' @export
gene_smooth <- function(mat, bincells = 1, normalize = TRUE, ncore = 1){
  # Kalman filter on DLM model
  mat <- as.matrix(mat)
  if (bincells > 1){
    mat <- binExp(mat, binsize = bincells)
  }
  dname <- dimnames(mat)
  if (ncore == 1){
    sm <- t(apply(mat, 1, dlm_helper))
  }
  else{
    message("Using ", ncore, " cores.")
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    sm <- foreach(i=1:nrow(mat), .packages=c("dlm"), .export=c("dlm_helper"), .combine=rbind) %dopar% {
      dlm_helper(mat[i, ])
    }
  }
  if (normalize){
    sm <- t(apply(sm, 1, function(g){
      mmin <- min(g)
      mmax <- max(g)
      (g - mmin)/(mmax - mmin)
    }))
  }
  dimnames(sm) <- dname
  
  return(as.data.frame(sm))
}

#' Compare gene expression
#' 
#' Compare gene expression.
#' @importFrom ggpubr rremove 
#' @import SingleCellExperiment
#' @import RColorBrewer
#' @import ggplot2
#' @param sce Expression matrix.
#' @param genes Use the mean of bincells for smoothing.
#' @param trunk Normalize the smoothed value to 0-1.
#' @param ncore ncores for parallel.
#' @return Smoothing expression.
#' @export
compare_expr <- function(sce, 
                           genes,
                           trunk = NULL, 
                           branch = NULL,
                           color = NULL,
                           linestyle = NULL,
                           datatype = "logcounts",
                           stateid = "branch",
                           trajectory = "trajectory",
                           lw = 2,
                           span = 0.6,
                           ...){
  anno <- colData(sce)
  sce <- sce[genes, anno[,stateid] %in% c(trunk, unlist(branch))]
  brs <- list()
  for (br in names(branch)){
    st <- c(trunk, unlist(branch[[br]]))
    mask <- colData(sce)[,stateid] %in% st
    sub_sce <- sce[, colData(sce)[,stateid] %in% st]
    sub_sce <- sub_sce[, order(colData(sub_sce)[,trajectory])]
    
    clnum <- table(colData(sub_sce)[, stateid])
    sub_sce$position <- c(seq(0, 0.5, length.out = sum(clnum[trunk])), seq(0.5+1/ncol(sub_sce), 1, length.out = sum(clnum[unlist(branch[[br]])])))
    sub_sce$trend <- br
    
    sm <- gene_smooth(assay(sub_sce, datatype), normalize = FALSE, ...)
    # TODO global application of bin expression 
    # sm[apply(sm, 1, max) < 2 | rowSds(as.matrix(sm))/rowMeans(sm) < 0.1, ] <- 0
    brs[[br]] <- as.data.frame(cbind(colData(sub_sce), t(sm)))
    
    ctexp <- do.call(rbind, brs)
    
    # Normalize
    ctexp[,genes] <- apply(ctexp[,genes], 2, function(g){
      (g - min(g))/(max(g) - min(g))
    })
  }
  
  # Ploting
  plist <- list()
  for (g in genes){
    ctexp$Expression <- ctexp[, g]
    ctexp$trend <- factor(ctexp$trend, levels = names(branch))
    plist[[g]] <- ggplot(ctexp, aes(x = position, y = Expression, colour = trend)) +
      geom_vline(xintercept = 0.5, linetype = "dashed") + 
      stat_smooth( span = span, se = FALSE, size = lw, aes(linetype=trend)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()) +
      labs(x = NULL, y = NULL, color = "Branch", title = g) + 
      rremove("legend")
    if (!is.null(color)){
      plist[[g]] <- plist[[g]] + scale_color_manual(values = color)
    }
    if (!is.null(linestyle)){
      plist[[g]] <- plist[[g]] + scale_linetype_manual(values=linestyle)
    }
  }
  ggarrange(plotlist = plist)
}

#' Penalty correlation
#' 
#' Penalty correlation
#' @importFrom stats cor
#' @param mat Gene expression matrix, rows are genes and cols are cells.
#' @param time Pseudotime
#' @param changepoint The threshold to define change point of a gene in pseudotime
#' @param beta Penalty coeffeicent
#' @return Corrected correlation matrix and change point.
#' @export
pcor = function(mat, 
                time, 
                changepoint=0.5,
                beta=1){
  # Normalization
  # mat = t(apply(mat, 1, function(g){
  #   mmin = min(g)
  #   mmax = max(g)
  #   (g - mmin)/(mmax - mmin)
  # }))
  pseudotime = (time-min(time))/(max(time)-min(time))
  
  # Identify change time point
  o = order(pseudotime)
  mat = mat[, o]
  pseudotime = pseudotime[o]
  e = apply(mat, 1, function(g){
    point = c()
    for (t in 1:length(pseudotime)){
      if (s(g, pseudotime) > 0){
        if (g[t] > changepoint){
          point = c(point, pseudotime[t])
          break
        }
      }
      else{
        if (g[t] < changepoint){
          point = c(point, pseudotime[t])
          break
        }
      }
      
    };point
  })
  
  # Calculate correlation
  pr = cor(t(mat))
  
  # Coefficient
  pcor_mat = matrix(rep(0, dim(pr)[1]^2), nrow = dim(pr)[1])
  dimnames(pcor_mat) = dimnames(pr)
  for (i in 1:nrow(pr)){
    for (j in 1:ncol(pr)){
      if (e[i] < e[j]){
        pcor_mat[i, j] = f(mat[i,], mat[j,], pr[i,j], beta)*pr[i,j]
      }
    }
  }
  return(list(pcor_mat = pcor_mat, changepoint = e))
}

s = function(g, pseudotime){
  return(sign(cor(as.numeric(g), as.numeric(pseudotime))))
}

f = function(x1, x2, pr, beta){
  n = length(x1)
  if (pr >= 0){
    return((1 - sum(abs(x1-x2))/n)^beta)
  }
  else{
    return((1 - sum(abs(x1+x2-1))/n)^beta)
  }
}

binExp = function(mat, binsize = 10){
  cells = colnames(mat)
  cellnum = length(cells)
  binum = cellnum %/% binsize
  rd = cellnum - binum * binsize
  set.seed(1)
  mat = mat[, !cells %in% sample(cells, rd)]
  bin_id = rep(paste0("bin", 1:binum), each = binsize)
  bin_id = factor(bin_id, levels = paste0("bin", 1:binum))
  m = aggregate(t(mat), list(bin_id), mean)
  m = t(m[,-1])
  colnames(m) = paste0("bin", 1:ncol(m))
  return(m)
}