
#########################################################################
# ------------------------- bespoke_distance() --------------------------

# Author: Kate Pereira Maia
# Date: 06/2021

# bespoke_distance: Calculates distances in interaction patterns between 
# insect species based on two arguments: data and dist.
# data: Considers three types of interaction data for distance calculations:
# binary (1 and 0), raw (interaction frequecies) or sumto1 (with abundances
# removed - see remove abundance function).
# dist: Calculates either, Euclidian, Gower, Bray-Curtis or Jaccard distances.

# ATTENTION: NOT ALL COMBINATIONS MAKE SENSE

bespoke_distance <- function(web, data = c("bin", "raw", "sumto1"), dist = c("eucl", "gow", "bray", "jac")) {
  
  if (data == "bin") {web[web != 0] <- 1}
  if (data == "sumto1") {web <- apply(web, 2, remove_abundance)}
  
  # diag = T | F are equal 0
  if (dist == "bray") { 
    web <- vegdist(t(web), method = "bray", diag = T, upper = F)
  } else if (dist == "jac") {
    web <- vegdist(t(web), method = "jaccard", diag = T, upper = F)
  } else if (dist == "euc") {
    web <- vegdist(t(web), method = "euclidean", diag = T, upper = F)
  } else {web <- vegdist(t(web), method = "gower", diag = T, upper = F)}
    
  return(web)
}
