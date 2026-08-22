# Small internal graph helpers (vendored from the paper code). Not exported.

decimal_to_dag <- function(decimal, n_nodes) {
  binary <- as.integer(intToBits(decimal))[1:(n_nodes * (n_nodes - 1))]
  matrix <- matrix(0, n_nodes, n_nodes)
  idx <- 1
  for (i in 1:n_nodes) {
    for (j in 1:n_nodes) {
      if (i != j) {
        matrix[i, j] <- binary[idx]
        idx <- idx + 1
      }
    }
  }
  return(matrix)
}

# Check that a 0/1 adjacency matrix is acyclic.
is_dag <- function(adj_matrix) {
  n <- nrow(adj_matrix)
  visited <- logical(n)
  rec_stack <- logical(n)
  dfs <- function(v) {
    visited[v] <<- TRUE
    rec_stack[v] <<- TRUE
    for (i in 1:n) {
      if (adj_matrix[v, i] == 1) {
        if (!visited[i]) {
          if (dfs(i)) return(TRUE)
        } else if (rec_stack[i]) {
          return(TRUE)
        }
      }
    }
    rec_stack[v] <<- FALSE
    return(FALSE)
  }
  for (i in 1:n) {
    if (!visited[i]) {
      if (dfs(i)) return(FALSE)
    }
  }
  return(TRUE)
}

# Check that each node has at most one parent (a directed tree/forest).
max_one_parent <- function(adj_matrix) {
  return(all(colSums(adj_matrix) <= 1))
}
