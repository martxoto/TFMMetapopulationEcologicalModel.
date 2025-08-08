# robustness_rw_ls() returns the robustness of the network, and also those of the different guilds
# plants do not get extinct because they are not predated or pollinated anymore.
# insects get extinct if one of their life stages lose all its interactions.
# A species / life stage is considered to be extinct if more than THRS_EXT interaction events have disapeared.
# Each time a species has its interactions rewired, a proportion THRS_FLEX get redistributed to other species in proportion to the abundance of the alternative resources/preys.
# additionally, it returns the proportion of remaining species against the proportion of extinct plant species so as to plot graphs similar to Figure 1.
# extinction scenario: preferentially from the most common to the least common plant species

robustness_rw_ls_weighted <- function(AObs, SP, BSP_C, LS, THRS_EXT, THRS_FLEX, scenario){
  # AObs is the 2D/3D array describing adjacencies among species
  # SP is a matrix object, which rows are named after species, and columns with guilds. SP[i, j] = 1 if species i belongs to guild j
  # BSP_C is a vector of plant commonness in the landscape
  # THRS_EXT is the proportion of interaction events that can disapear before a species is considered extinct
  # THRS_FLEX is the proportion of interaction events that can be reallocated after primary extinction
  # scenario is a character string = {"random", "mostcommon", "leastcommon"}; "random" stands for a random removal scenario while "mostcommon" refers to the most common species to the least scenario, and "leastcommon" the other way around
  # LS is a 3D array describing the life stage of each species when involved in a given type of interactions
  # it has the same dimension than SP, but elements are = "larva" or "adult" for insects
  
  if (is.array(AObs) == FALSE){
    stop("AObs must be a 2D or 3D array.")
  }
  if (dim(AObs)[1] != dim(AObs)[2]){
    stop("Each AObs[,, k] must be square matrices.")
  }
  if (length(dim(SP)) != 2){
    stop("SP must be a matrix, which rows are named after species, and columns with guilds.")
  }
  if (nrow(SP) != dim(AObs)[1]){
    stop("SP must have as many rows as AObs[,, k].")
  }
  if (!is.array(LS)){
    stop("LS must be a 3D array.")
  }
  if (length(dim(LS)) != 3){
    stop("LS must be a 3D array.")
  }
  else {
    if (length(intersect(colnames(LS[1,,]), c("larva", "adult"))) != 2){
      stop("LS must describe insects life stage in its third dimension.")
    }
  }
  if (((!is.numeric(THRS_EXT)) | (THRS_EXT > 1) | (THRS_EXT < 0)) | ((!is.numeric(THRS_FLEX)) | (THRS_FLEX > 1) | (THRS_FLEX < 0))){
    stop("THRS_EXT and THRS_FLEX must be numeric, positive and below 1.")
  }
  if (is.character(scenario) == FALSE){
    stop("scenario must be a character string. It is equal to 'random' or 'mostcommon'.")
  }
  if ((scenario != "random") & (scenario != "mostcommon") & (scenario != "leastcommon")){
    stop("scenario must be a character string. It is equal to 'random', 'mostcommon', or 'leastcommon'.")
  }
  
  N_tot <- nrow(SP) # total number of species
  N_pl <- sum(SP[, 'PL']) # number of plant species
  Species <- rownames(SP)
  pl_list <- Species[which(SP[, 'PL'] != 0)]
  
  R <- c() # robustness of all guilds
  
  scd_extinctions <- matrix(0, N_tot, 0); rownames(scd_extinctions) <- rownames(SP)
  rm_seq <- c()
  ins_rm_seq <- c()
  
  if (scenario == "mostcommon"){
    p_rm <- BSP_C/sum(BSP_C) # probability to be removed based on species commonness (from most to least common)
  }
  if (scenario == "leastcommon"){
    p_rm <- (1/BSP_C); p_rm <- p_rm/sum(p_rm) # probability to be removed based on species commonness (from least to most common)
  }
  if (scenario == "random") {
    p_rm <- rep(1, N_pl); p_rm <- p_rm/N_pl
    names(p_rm) <- pl_list
  }
  
  ARef <- AObs # save the observed matrix
  # amount of interaction events for each insect and at its two life stages
  ResRefAdult <- apply(ARef, 3, colSums)[(N_pl+1):N_tot,]*LS[,, "adult"]; ResRefAdult <- rowSums(ResRefAdult)
  ResRefLarva <- apply(ARef, 3, colSums)[(N_pl+1):N_tot,]*LS[,, "larva"]; ResRefLarva <- rowSums(ResRefLarva)
  ResRef <- cbind(ResRefLarva, ResRefAdult); colnames(ResRef) <- c("larva", "adult")

  # extinction sequence
  while (length(rm_seq) != N_pl){ # while there are plants to remove
    p_rm[rm_seq] <- 0; p_rm <- p_rm/sum(p_rm) # remove already removed species
    PL_RM <- my_sample(pl_list, 1, repl = FALSE, pi = p_rm) # pick a plant species to remove based on its commonness
    rm_seq <- c(rm_seq, PL_RM) # record the species name to remove
    rm_links <- which(AObs[PL_RM, , ] != 0, arr.ind = TRUE) # which interactions could be lost because of plant extinction? (consumer/predator as 'row') in which network type ? (as 'col')
    if (THRS_FLEX != 0){
      int_with_pl_rm <- AObs[PL_RM, , ]
      rm_links <- cbind(rm_links, int_with_pl_rm[int_with_pl_rm != 0])
      
      ABin <- AObs; ABin[ABin != 0] <- 1
      RO_g <- res_overlap_graph(ABin) # builds resource-overlap graph
      RW_g <- rewiring_graph(ABin, RO_g) # builds rewiring graph
      
      AObs <- rewire_weighted(rm_links, PL_RM, RW_g, AObs, LS, ResRef, THRS_FLEX, THRS_EXT)
    }
    else {
      AObs[PL_RM, , ] <- 0
    }

    ExtinctSp <- c()
    if (ncol(scd_extinctions) != 0){
      ExtinctSp <- which(rowSums(scd_extinctions) == 1)
    }
    Ins2Rm <- whichSpExtinct(AObs, THRS_EXT, ResRef, LS, ExtinctSp)
    AObs[Ins2Rm, , ] <- 0; AObs[, Ins2Rm, ] <- 0
    new_extinct <- vector('numeric', nrow(scd_extinctions))
    new_extinct[Ins2Rm] <- 1; new_extinct[which(Species == PL_RM)] <- 1
    scd_extinctions <- cbind(scd_extinctions, new_extinct)
  }

  extinct <- colSums(scd_extinctions)
  extinct <- cumsum(extinct)
  p_extinct <- extinct/N_tot
  p_remain <- 1-p_extinct; p_remain <- c(1, p_remain)
  for (i1 in colnames(SP)){
    sp_i1 <- which(SP[, i1] == 1) # species of type i1
    spi1_extinct <- scd_extinctions[sp_i1, ]
    if (length(sp_i1) > 1){
      spi1_extinct <- colSums(scd_extinctions[sp_i1, ])
    }
    spi1_extinct <- cumsum(spi1_extinct)
    pi1_extinct <- spi1_extinct/(sum(SP[, i1]))
    pi1_remain <- 1-pi1_extinct; pi1_remain <- c(1, pi1_remain)
    p_remain <- rbind(p_remain, pi1_remain)
  }
  rownames(p_remain) <- c('All species', colnames(SP))

  pl_removed <- c(0:N_pl); pl_removed <- pl_removed/N_pl
  R <- vector('numeric', nrow(p_remain)); names(R) <- rownames(p_remain); R['PL'] <- NA

  c <- 1
  for (i1 in rownames(p_remain)){
    if (i1 != 'PL'){
       ext_curve <- splinefun(pl_removed, p_remain[i1, ]) # interpolate the extinction 'curve'
       ext_area <- integrate(ext_curve, 0, 1) # calculate the area under the curve = robustness
       R[i1] <- as.numeric(ext_area[[1]])
       c <- c+1
    }
  }
  
  NAGuild <- colnames(SP)[colSums(SP)==0]
  R[NAGuild] <- NA; p_remain[NAGuild, ] <- NA
  output_R <- list(robustness = R, remaining_sp = p_remain)
  return(output_R)
}