#' Accurate tiling of spatial single-cell data
#'
#' An algorithm for segmenting single-cell resolution
#' spatial omics data into small multicellular tiles whose edges track with
#' natural tissue boundaries. These tiles can then be used in downstream
#' analysis to label and define tissue regions across samples.
#'
#' @name tessera
#' @useDynLib tessera
#' @import Rcpp
#' @importClassesFrom Matrix dgCMatrix 
#' @importFrom methods as is
#' @importFrom Matrix Matrix sparseMatrix
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom utils head
#' @importFrom Rcpp evalCpp sourceCpp loadModule
#' @importFrom rlang .data
#' @importFrom sf st_cast
#' @importFrom furrr future_map
"_PACKAGE"
