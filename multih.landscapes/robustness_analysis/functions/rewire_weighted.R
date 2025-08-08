# rewire_weighted() returns the graph A (2D or 3D array) after rewiring between two trophic levels when possible
# A rewiring is done based on the rewiring graph (RW).
rewire_weighted <- function(int2rw, sp2rm, RW, A, InsLS, ObsRes, FLEX, EXT){
  # int2rw is the list of interactions to rewire or remove if can't be rewired (a matrix with 3 columns), and their weigth
  # 1st column = predators/consumers of the removed species; 2nd column = interaction type ID
  # note that int2rw lists only the interactions involving the lower taxon deleted by the robustness analysis
  
  # RW is the rewiring graph (2D or 3D array)
  # A is the 2D/3D array describing adjacencies among species
  # InsLS is a 3D array describing the life stage of each species when involved in a given type of interactions
  # ObsRes is a matrix describing the number of interaction events a species has with its resources/preys/hosts for each life stage
  # FLEX is the proportion of interaction to reallocate
  # EXT is the proportion of interaction events that can disapear before a species is considered extinct
  
  if (is.matrix(int2rw) == FALSE){
    stop("int2rw must be a 3-columns matrix.")
  }
  else {
    if (ncol(int2rw) != 3){
      stop("int2rw must be a 3-columns matrix.")
    }
  }
  
  if ((is.array(A) == FALSE) | (is.array(RW) == FALSE)){
    stop("A and RW must be arrays.")
  }
  if (sum(abs(dim(A)-dim(RW))) != 0){
    stop("A and RW must have the same dimensions.")
  }
  if (((!is.numeric(EXT)) | (EXT > 1) | (EXT < 0)) | ((!is.numeric(FLEX)) | (FLEX > 1) | (FLEX < 0))){
    stop("EXT and FLEX must be numeric, positive and below 1.")
  }
  if (!is.array(InsLS)){
    stop("InsLS must be a 3D array.")
  }
  if (length(dim(InsLS)) != 3){
    stop("InsLS must be a 3D array.")
  }
  else {
    if (length(intersect(colnames(InsLS[1,,]), c("larva", "adult"))) != 2){
      stop("InsLS must describe insects life stage in its third dimension.")
    }
  }
  if (!is.matrix(ObsRes)){
    stop("ObsRes must be a matrix which rows are species in the network, and columns life stages.")
  }
  if (length(intersect(colnames(ObsRes), c("larva", "adult"))) != 2){
    stop("ObsRes must have two columns named 'larva' and 'adult'.")
  }
  
  
  Arw <- A # create and initialise the adjacencies with rewiring
  
  if (nrow(ObsRes) != nrow(InsLS)){
    stop("ObsRes must be a matrix which rows are species in the network.")
  }
  SPl <- nrow(A)-nrow(ObsRes) # number of plants
  STot <- nrow(A) # total number of species
  
  int_WT <- unique(int2rw[, 2]) # type of interactions to remove
  for (i1 in int_WT){ # for each type of interaction
    IntType_i1 <- which(int2rw[, 2] == i1)
    int2rw_i1 <- int2rw[IntType_i1, ] # interactions to remove of type i1 (more precisely, predators/consumers to which a new resource should be attributed)
    for (i2 in 1:length(IntType_i1)) { # for each interaction of a given type
      if (length(IntType_i1) != 1){
        ut <- int2rw_i1[i2, 1] # upper taxon involved in interaction i2 of type i1
        NEvents <- int2rw_i1[i2, 3] # number of interction events which could be deleted
      }
      else {
        ut <- int2rw_i1[1]
        NEvents <- int2rw_i1[3]
      }
      overlap_res <- which(RW[, ut, i1] != 0) # can ut rewire its lost interaction? (list of consumers/predators with which ut share resources)

       if (length(overlap_res) != 0){ # if the answer is yes
         alternative <- Arw[, overlap_res, i1] # what are the alternative resources in network type i1?

         if (is.null(dim(alternative)) == FALSE){
           alternative <- rowSums(alternative)
         }
         alternative[sp2rm] <- 0 # remove sp2rm from this list of resources
        
         if (sum(alternative) != 0){ # if there's actual alternative resources/hosts
           new_res <- my_sample(alternative, NEvents, repl = TRUE, alternative/sum(alternative)) # pick new resources (as many as there are interaction events which could be deleted), with probability proportional to the number of other consumers/predators feeding on it
           new_res <- names(new_res) # get the name
           
           prop_res <- my_sample(c(TRUE, FALSE), length(new_res), repl = TRUE, c(FLEX, 1-FLEX))
           new_res <- new_res[prop_res]
           if (length(Arw[new_res, ut, i1]) != length(new_res)){
             stop()
           }
           Arw[new_res, ut, i1] <- Arw[new_res, ut, i1]+1 # create the alternative interaction
         }
        }

        # remove the deleted interaction
        Arw[sp2rm, ut, i1] <- 0
        Arw[ut, sp2rm, i1] <- 0

        # does species ut have enough resources to persist? (should they be deleted from the web?)
        InsResAdult <- apply(Arw, 3, colSums)[(SPl+1):STot,]*InsLS[,, "adult"]; InsResAdult <- rowSums(InsResAdult)
        InsResLarva <- apply(Arw, 3, colSums)[(SPl+1):STot,]*InsLS[,, "larva"]; InsResLarva <- rowSums(InsResLarva)

        if ((InsResAdult[ut-SPl] < EXT*ObsRes[ut-SPl, "adult"]) | (InsResLarva[ut-SPl] < EXT*ObsRes[ut-SPl, "larva"])){ # if ut should be removed
          int2rw_ut <- which(Arw[ut, , ] != 0, arr.ind = TRUE) # outward links for species ut
          int_with_ut <- Arw[ut, , ]
          if (nrow(int2rw_ut) != length(int_with_ut[int_with_ut != 0])){
            print(int2rw_ut)
            print(int_with_ut[int_with_ut != 0])
            stop()
          }
          int2rw_ut <- cbind(int2rw_ut, int_with_ut[int_with_ut != 0])
          Arw <- rewire_weighted(int2rw_ut, ut, RW, Arw, InsLS, ObsRes, FLEX, EXT)
        }
     }
  }
  
  return(Arw)
  
}