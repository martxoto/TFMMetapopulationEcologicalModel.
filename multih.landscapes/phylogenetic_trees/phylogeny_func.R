
##################################################################################################
#------------------------------ Function PD: Phylogenetic diversity ------------------------------
# Author: Kate Pereira Maia
# Date: 19/04/2020

# Task: Constructs the phylogenetic tree and calculates the phylogenetic diversity - measured as the sum of branch lengths - of a plant community. The tree is a chop of the Daphne tree. 

# Arguments:
# i. plant_comm: the vector of plant species to be included in the tree
# ii. solution: a dataframe listing problems/solutions for species not available in the phylogeny (Daphne)
# iii. tree: the large phylogenectic tree (Daphne) to be chopped

# Returns: A list with 3 elments:
# i. A dataframe with 3 (or 4, see obs) values:
#  - n_tree: number of alternative trees - as some species require multiple solutions - which were built for the plant_commm 
#  - m_PD: mean phylogenetic diversity of the alternative trees built
#  - sd_PD: standard deviation of the mean phylogenetic diversity of the alternative trees built
#  - obs: Whether ranunculus_acris (typo) is part of the community - later fixed in the dataset
# ii. A vector with values of PF for the alternative trees
# iii. The community dataframe (comm_df): the tree positions used to build the alternative trees

# GO TO THE END (LINE 185) FOR AN EXAMPLE OF HOW TO USE THE FUNCTION PD

#----------------------------- Uploading packages ----------------------------

library(ape); library(tidyverse)

#-------------------- Reading data and creating arguments --------------------

interac_data <- read.table("./Data/all_web_interactions_final.txt", header = TRUE, sep = "\t")
interac_data <- interac_data[!interac_data$Web %in% c("CP_P", "LM_P"), ] # only plant interactions
interac_data$Site <- factor(interac_data$Site); levels(interac_data$Site)[21] <- "Pembrey_Beach"
interac_data$Lower_Taxon <- factor(interac_data$Lower_Taxon); levels(interac_data$Lower_Taxon)[1] <- "moss"
plant_comm <- levels(interac_data$Lower_Taxon)

tree <- read.tree("./Data/DaPhnE_01.tre")

solution <- read.table("plant_solution.txt", header = TRUE, sep = "\t", na.strings = "")

#-------------------- Reading data and creating arguments --------------------

PD <- function(plant_comm, solution, tree) {
  
  # Community dataframe to be filled with the tree position of all species in plant_comm
  comm_df <- data.frame(sp = plant_comm, S1 = NA)
  
  # Deals with Ranunculus_acris/ranunculus_acris
  if (any("ranunculus_acris" == comm_df$sp)) {racris <- TRUE} else {racris <- FALSE}
  if (any("Ranunculus_acris" == comm_df$sp)) {Racris <- TRUE} else {Racris <- FALSE}
  
  #----------------------- OK plants ----------------------
  
  # ok_comm is the ok subset of plant_comm
  # used to create ok_df: col1 = ok_comm, col2 = position of ok_comm in Daphne
  # position in ok_df is then matched to comm_df
  # THIS STRUCTURE WILL REPEAT ITSELF FOR PROBLEM PLANTS
  
  ok_comm <- plant_comm[!plant_comm %in% solution$miss_plant_list]
  ok_df <- data.frame(sp = grep(paste(ok_comm, collapse = "|"), tree$tip.label, value = TRUE), 
                      posit = grep(paste(ok_comm, collapse = "|"), tree$tip.label))
  
  if (length(ok_comm) != nrow(ok_df)) {print("Lenght OK")}
  
  comm_df$S1 <- ok_df[match(comm_df$sp, ok_df$sp), "posit"]
  
  prob_table <- table(solution$solution) == 0
  
  #-------------------- PROBLEM plants --------------------
  
  solution <- solution[solution$miss_plant_list %in% plant_comm,]
  solution <- cbind(gn = unlist(lapply(strsplit(as.character(solution$miss_plant_list), "_"), 
                                       function(x) x <- x[1])), solution)
  
  # there are potentially 4 types of problems to be sorted
  table(solution$solution)
  
  # alternative sp with one alternative: can be dealt with as typo sp
  alt_to_typo <- table(solution[solution$solution == "alternative", "gn"])
  alt_to_typo <- alt_to_typo[alt_to_typo == 1]
  solution$solution[match(names(alt_to_typo), solution$gn)] <- "typo"
  
  
  prob <- table(solution$solution) > 0
  prob_table[match(names(prob), names(prob_table))] <- prob
  prob <- as.data.frame(t(prob_table))
  
  #----- 1) TYPO: match to a different name in tree -----
  
  if (prob$typo) {
    typo_comm <- as.character(solution[solution$solution == "typo", "which"])
    typo_df <- data.frame(sp = grep(paste(typo_comm, collapse = "|"), tree$tip.label, value = TRUE), 
                          posit = grep(paste(typo_comm, collapse = "|"), tree$tip.label))
    
    if (length(typo_comm) != nrow(typo_df)) {print("Length typo")}
    
    typo_comm[!typo_comm %in% typo_df$sp]
    typo_df$sp[!typo_df$sp %in% typo_comm]
    
    # matching typo_df to typo_comm is equivalent to matching it to solution
    typo_df <- typo_df[match(typo_comm, typo_df$sp),]
    all(typo_comm == typo_df$sp)
    
    typo_comm <- as.character(solution[solution$solution == "typo", "miss_plant_list"])
    where <- which(!is.na(typo_df[match(comm_df$sp, typo_comm), "posit"]))
    what <- na.omit(typo_df[match(comm_df$sp, typo_comm), "posit"])
    comm_df$S1[where] <- what
  }
  
  # 2)----- GENUS: sample one from available species -----
  
  if (prob$genus) {
    genus_comm <- as.character(solution[solution$solution == "genus", "which"])
    genus_df <- data.frame(sp = grep(paste(genus_comm, collapse = "|"), tree$tip.label, value = TRUE), 
                           posit = grep(paste(genus_comm, collapse = "|"), tree$tip.label))
    genus_df <- cbind(genus = unlist(lapply(strsplit(as.character(genus_df$sp), "_"), function(x) x <- x[1])), 
                      genus_df)
    genus_df <- tapply(genus_df$posit, genus_df$genus, sample, 1)
    
    names(genus_df) <- solution[match(names(genus_df), solution$which), "miss_plant_list"]
    comm_df[match(names(genus_df), comm_df$sp), "S1"] <- genus_df
  }
  
  # 3)----- ALTERNATIVE: creates all scenarios -----
  
  if (prob$alternative) {
    tab_gn <- table(factor(solution[solution$solution == "alternative", "gn"]))
    n_scen <- prod(tab_gn)
    
    alt_comm <- as.character(solution[solution$solution == "alternative", "which"])
    alt_df <- data.frame(sp = grep(paste(alt_comm, collapse = "|"), tree$tip.label, value = TRUE), 
                         posit = grep(paste(alt_comm, collapse = "|"), tree$tip.label))
    alt_df <- cbind(gn = unlist(lapply(strsplit(as.character(alt_df$sp), "_"), 
                                       function(x) x <- x[1])), alt_df)
    
    if (length(tab_gn) == 1) {scen <- expand.grid(1:3)}
    if (length(tab_gn) == 2) {scen <- expand.grid(1:3, 1:3)}
    if (length(tab_gn) == 3) {scen <- expand.grid(1:3, 1:3, 1:3)}
    colnames(scen) <- names(tab_gn)
    for (i in 1:ncol(scen)) {
      gn <- colnames(scen)[i]
      sp <- alt_df[alt_df$gn == gn, "posit"]
      scen[,i] <- sp[scen[,i]]
    }
    
    scen_comm_df <- comm_df[,rep(2, each = (n_scen - 1))]
    colnames(scen_comm_df) <- paste("S", 2:n_scen, sep = "")
    comm_df <- cbind(comm_df, scen_comm_df)
    
    colnames(scen) <- unique(solution$miss_plant_list[solution$solution == "alternative"])
    scen <- t(scen)
    comm_df[match(rownames(scen), comm_df$sp), -1] <- scen
  } else {n_scen <- 1}
  
  # 4)----- N: generally Poaceae (in PL_CP e PL_LM) - CURRENTLY REMOVING -----
  
  if (prob$N) {
    comm_df <- comm_df[-which(is.na(comm_df$S1)),]
    if (sum(is.na(comm_df)) != 0) {print("NA")}
  }
  
  #-------------------- Construction and calculations --------------------
  
  if (Racris & racris) {comm_df <- comm_df[-which(comm_df$sp == "ranunculus_acris"),]}
  
  phyl_div <- c()
  for (i in 2:(n_scen + 1)) {
    to.keep <- comm_df[,i]
    to.drop <- seq(1:length(tree$tip.label))[!seq(1:length(tree$tip.label)) %in% to.keep]
    comm_tree <- drop.tip(tree, to.drop)
    if (length(comm_tree$tip.label) != length(to.keep)) {print("PHYL")}
    pd <- sum(comm_tree$edge.length)
    phyl_div <- c(phyl_div, pd)
  }
  
  ret_summ <- data.frame(n_tree = length(phyl_div), mean_PD = mean(phyl_div), sd_PD = sd(phyl_div))
  
  if (racris) {ret_summ <- cbind(ret_summ, r.acris = TRUE)}
  
  to_return <- list(ret_summ, phyl_div, comm_df)
  return(to_return)
  
}

#--------------------------- Example of PD Function --------------------------

for (i in 1:300) {
  print(i)
  test <- plant_comm[sample(1:length(plant_comm), 150)]
  result <- PD(plant_comm = test, solution, tree)
}

#--------------------- Running for empirical communities ---------------------

site_data <- read.table("./Data/site_level_data.csv", header = TRUE, sep = ",")
site_data$phyl_div <- NA
site_data$Site <- gsub(" ", "_", site_data$Site)
site <- site_data$Site

for (i in 1:length(site)) {
  plant_comm_site <- unique(interac_data[interac_data$Site == site[i], "Lower_Taxon"])
  result <- PD(plant_comm = plant_comm_site, solution, tree)
  site_data[i, "phyl_div"] <- result[[1]][2]
}

site_data %>% 
  mutate(MDT = recode(MDT, "1" = "Monad", "2" = "Dyad", "3" = "Triad")) %>%
  mutate(MDT = factor(MDT, levels = c("Monad", "Dyad", "Triad"))) %>% 
  ggplot(aes(x = MDT, y = phyl_div)) +
  geom_boxplot() +
  xlab("") + ylab("Phylogenetic diversity") +
  theme_classic()

#write.table(site_data, "./Data/site_level_data.csv", sep = ",")
