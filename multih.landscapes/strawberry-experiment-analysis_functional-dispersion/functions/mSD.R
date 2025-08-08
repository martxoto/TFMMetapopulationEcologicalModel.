
#########################################################################
# -------------------------------- mSD() --------------------------------

# Author: Kate Pereira Maia
# Date: 06/2021

# mSD: Code from Maire et al which calculates the mean squared deviation 
# between the original distances (provided by the argument original_dist) 
# and the euclidian distances in the multidimensional space created with 
# PCoA. It returns a vector with the mSD values for multidimensional spaces
# with 2 up to nbdim (I guess it is the no of orthogonal PCoA axes) axes.

mSD <- function(original_dist){
  
  dist_raw <- dist_st <- list()
  
  pcoa <- pcoa(as.dist(original_dist), correction = "cailliez")
  nbdim <- ncol(pcoa$vectors)
  coord <- pcoa$vectors[, 1:nbdim]
  colnames(coord) <- paste("PC", 1:nbdim, sep = "")
  
  for (k in 2:nbdim) {
    # Calculates the euclidian distance up to the kth PC and saves in dist_kD
    eval(parse(text = paste("dist_", k ,"D <- dist(coord[, 1:", k, "], method = 'euclidean')", sep = "")))
    
    # Saves dist_kD as an element of dist_raw under the name m_kD
    eval(parse(text = paste("dist_raw$m_", k, "D <- dist_", k, "D", sep = "")))
  }
  
  meanSD <- seq(nbdim - 1)
  names(meanSD) <- c(paste("m_", 2:nbdim, "D", sep = ""))
  S <- ncol(as.matrix(original_dist)) # number of flower visitors
  x <- original_dist
  
  for (k in 2:nbdim) {
    # y = dist_kD (see above): euclidian distance matrix for the kth dimension
    eval(parse(text = paste("y <- dist_", k, "D", sep = "")))
    
    # Puts y (functional space) in the same scale as x (original distances)
    yst <- y / max(y) * max(x) 
    
    # Saves dist_kD as an element of dist_st under the name m_kD
    eval(parse(text = paste("dist_st$m_", k, "D<-dist_", k, "D", sep = ""))) # WHY?
    
    meanSD[paste("m_", k, "D", sep = "")] <- round(((sum((x - yst) ^ 2))/(S * (S - 1) / 2)), 6)
  }
  return(meanSD)
}