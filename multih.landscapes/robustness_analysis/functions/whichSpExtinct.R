# whichSpExtinct() returns a vector of species ID which are considered extinct according to the following rules:
# A species can be considered extinct because it doesn't have enough resources to feed on anymore.
# The function needs a threshold for extinction as one of its input argument. This threshold corresponds to
# the proportion of interaction events below which the species is considered not to have enough resources.

whichSpExtinct <- function(Web, ExtThrs, ObsInt, LifeStages, ExtinctSp){
  # Web is the 2D/3D array describing adjacencies among species
  # ExtThrs is the proportion of interaction events that can disapear before a species is considered extinct
  # ObsInt is a matrix describing the number of interaction events a species has with its resources/preys/hosts for each life stage
  # LifeStages is a 3D array describing the life stage of each species when involved in a given type of interactions
  # ExtinctSp is a vector of species that are already extinct
  
  if (!is.array(Web)){
    stop("Web must be a 2D or 3D array.")
  }
  if (dim(Web)[1] != dim(Web)[2]){
    stop("Each Web[,, k] must be square matrices.")
  }
  if ((!is.numeric(ExtThrs)) | (ExtThrs > 1) | (ExtThrs < 0)){
    stop("ExtThrs must be numeric, positive and below 1.")
  }
  if (!is.matrix(ObsInt)){
    stop("ObsInt must be a matrix which rows are species in the network, and columns life stages.")
  }
  if (length(intersect(colnames(ObsInt), c("larva", "adult"))) != 2){
    stop("ObsInt must have two columns named 'larva' and 'adult'.")
  }
  if (!is.array(LifeStages)){
    stop("LifeStages must be a 3D array.")
  }
  if (length(dim(LifeStages)) != 3){
    stop("LifeStages must be a 3D array.")
  }
  else {
    if (length(intersect(colnames(LifeStages[1,,]), c("larva", "adult"))) != 2){
      stop("LifeStages must describe insects life stage in its third dimension.")
    }
  }
  if ((length(ExtinctSp) != 0) & (!is.vector(ExtinctSp))){
    stop("ExtinctSp must be a vector.")
  }
  
  NSp <- nrow(Web) # number of species
  NPl <- NSp-nrow(ObsInt) # number of plants
  
  if (nrow(ObsInt) != nrow(LifeStages)){
    stop("ObsInt must be a matrix which rows are species in the network.")
  }


  # amount of interaction events for each insect and at its two life stages
  ResAsAdult <- apply(Web, 3, colSums)[(NPl+1):NSp,]*LifeStages[,, "adult"]; ResAsAdult <- rowSums(ResAsAdult)
  ResAsLarva <- apply(Web, 3, colSums)[(NPl+1):NSp,]*LifeStages[,, "larva"]; ResAsLarva <- rowSums(ResAsLarva)

  NewExtinctSp <- c()
  for (Sp in (NPl+1):NSp){
    TestPersist <- (length(intersect(Sp, ExtinctSp)) == 0) # is Sp already extinct?
    if (TestPersist){
      if ((ResAsAdult[Sp-NPl] < ExtThrs*ObsInt[Sp-NPl, "adult"]) | (ResAsLarva[Sp-NPl] < ExtThrs*ObsInt[Sp-NPl, "larva"])){
        NewExtinctSp <- c(NewExtinctSp, Sp)
      }
    }
  }

  return(NewExtinctSp)
}