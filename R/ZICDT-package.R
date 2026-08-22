#' ZICDT: Directed Tree Structure Learning for Zero-Inflated Compositional Nodes
#'
#' Each node is a whole composition on the probability simplex. The method scores
#' candidate edges with a KL divergence (finite on exact zeros), models the child
#' conditional mean as a mixture of a baseline and a parent-driven, column-stochastic
#' transition map, fits parameters by EM, and finds the globally optimal directed
#' tree/forest via Chu-Liu/Edmonds. The edge penalty is chosen by cross-validation.
#'
#' @keywords internal
#' @useDynLib ZICDT, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom igraph graph_from_adjacency_matrix V "V<-" layout_with_sugiyama
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply
#' @importFrom stats aggregate sd rgamma rmultinom runif
#' @importFrom graphics plot mtext
"_PACKAGE"
