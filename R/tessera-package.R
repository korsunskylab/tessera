#' tessera
#'
#' Spatial tiling 
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
