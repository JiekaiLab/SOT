% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudotime.R
\name{pcor}
\alias{pcor}
\title{Penalty correlation}
\usage{
pcor(mat, time, changepoint = 0.5, beta = 1)
}
\arguments{
\item{mat}{Gene expression matrix, rows are genes and cols are cells.}

\item{time}{Pseudotime}

\item{changepoint}{The threshold to define change point of a gene in pseudotime}

\item{beta}{Penalty coeffeicent}
}
\value{
Corrected correlation matrix and change point.
}
\description{
Penalty correlation
}
