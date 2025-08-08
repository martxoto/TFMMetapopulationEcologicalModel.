
#########################################################################
# ------------------------- remove_abundance() --------------------------

# Author: Kate Pereira Maia
# Date: 06/2021

# remove_abundance: Given the way that interaction data is collected, the
# interaction frequency of an insect species is equivalent to its abundance.
# This function removes insect abundances by turning absolute interaction 
# frequencies (colSums of interaction matrix) into relative interaction 
# frequencies (colSums of interaction matrix summing to 1). Therefore, 
# resulting matrix keeps interactions frequencies but removes insect 
# abundance.

remove_abundance <- function(x){sumx <- rep(sum(x), length(x)); newx <- x/sumx}

