# res_overlap_graph() returns the 2D/3D array describing which species share resources/preys in which type of network (slice)

res_overlap_graph <- function(G){
  # G is an (2D or 3D) array describing species interactions
  # if 2D, it is the adjacency matrix of the network, if 3D, each layer G[,, k] is an adjacency matrix
  
  if ((length(dim(G)) != 2) & (length(dim(G)) != 3)){
    stop("G must be a 2D or 3D array describing species interactions.")
  }
  S <- dim(G)[1] # number of species
  
  if (dim(G)[1] != dim(G)[2]){
    stop("Each slice of G must be a square matrix: the first two dimensions of G must be identical.")
  }
  
  RO_G <- array(0, dim = dim(G), dimnames = dimnames(G)) # creating and initialising the resource-overlap graph
  for (i1 in 1:S){
    degree_i1 <- colSums(G[i1,,]) # degree of species i1 in each layer (eg. network type)
    if (length(which(degree_i1 > 1)) != 0){ # if species i1 is shared by at least 2 predators/consumers in one of the network types
      shared_pred <- which(G[i1,,] != 0, arr.ind = TRUE) # list of predators/consumers sharing i1, second column tells in which situation
      slices <- unique(shared_pred[, 2]) # list of network slices where i1 is shared by at least 2 predators/consumers
      for (i2 in slices){ # for each slice
        shared_predi2 <- matrix(shared_pred[shared_pred[, 2] == i2, ], ncol = 2)# list of predators/consumers sharing i1 in slice i2
        RO_G[shared_predi2[, 1], shared_predi2[, 1], i2] <- RO_G[shared_predi2[, 1], shared_predi2[, 1], i2]+1
      }
    }
  }
  
  for (i3 in 1:dim(RO_G)[3]){
    diag(RO_G[,, i3]) <- 0
  }
  
  # At this stage, RO_G tells us whether two species share resources/preys, and how many
  RO_G[RO_G != 0] <- 1 # here, we only keep the information about whether resource sharing occurs or not.
  
  return(RO_G)
}