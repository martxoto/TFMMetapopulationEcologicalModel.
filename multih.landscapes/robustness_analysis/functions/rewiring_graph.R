# rewiring_graph() returns the 2D/3D array describing which species share resources/preys, and which feed on alternate ones
# RW_G[i, j] = 1 means that species i and j share at least one resource/prey, and that species i feeds on a resource/prey with which j doesn't interact
rewiring_graph <- function(G, RO_G){
  # G is an (2D or 3D) array describing species interactions
  # if 2D, it is the adjacency matrix of the network, if 3D, each layer G[,, k] is an adjacency matrix
  
  # RO_G is the resource overlap graph based on G
  
  if ((length(dim(G)) != 2) & (length(dim(G)) != 3)){
    stop("G must be a 2D or 3D array describing species interactions.")
  }
  S <- dim(G)[1] # number of species
  
  if (dim(G)[1] != dim(G)[2]){
    stop("Each slice of G must be a square matrix: the first two dimensions of G must be identical.")
  }
  
  if ((dim(RO_G)[1] != dim(G)[1]) | (dim(RO_G)[2] != dim(G)[2]) | (dim(RO_G)[3] != dim(G)[3])){
    stop("RO_G and G must have the same dimension.")
  }
  
  RW_G <- array(0, dim = dim(RO_G), dimnames = dimnames(RO_G)) # creating and initiating the rewiring graph
  slices <- seq(dim(RW_G)[3]) # the sequence of slices
  for (i1 in slices){
    overlap <- which((RO_G[,, i1] != 0) & (upper.tri(RO_G[,, i1])), arr.ind = TRUE)
    if (nrow(overlap) != 0){ # if some species share resources/preys
      for (i2 in 1:nrow(overlap)){
        sp1 <- overlap[i2, 1]; sp2 <- overlap[i2, 2]
        res_sp1 <- G[, sp1, i1]; res_sp2 <- G[, sp2, i1]
        comp_res <- res_sp1-res_sp2
        if (length(which(comp_res == 1)) != 0){
          RW_G[sp1, sp2, i1] <- 1
        }
        if (length(which(comp_res == -1)) != 0){
          RW_G[sp2, sp1, i1] <- 1
        }
      }
    }
  }
  
  return(RW_G)
}