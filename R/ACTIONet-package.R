#' ACTIONet: Robust multi-resolution analysis of single-cell datasets
#'
#' Joint archetypal/network-based analysis for multi-resolution decomposition of single-cell states.
#'
#' @section Usage: 
#' 
#' \enumerate{
#' \item ?reduce.sce to reduce raw count matrix in an SummarizedExperiment objects.
#' \item ?run.ACTIONet to run ACTIONet on SummarizedExperiment objects.
#' }
#' @section Useful links:
#' 
#' \enumerate{
#' \item Report bugs at \url{https://github.com/shmohammadi86/ACTIONet/issues}
#' \item Read the manuscript: 
#' \href{https://www.biorxiv.org/content/10.1101/746339v2}{A multiresolution framework to characterize single-cell state landscapes}.
#' }
#' 
#'
#' @name ACTIONet
#' @docType package
#' @author Shahin Mohammadi
#' @import Rcpp
#' @useDynLib ACTIONet, .registration=TRUE
#' @exportPattern ^[[:alpha:]]+ 
NULL
