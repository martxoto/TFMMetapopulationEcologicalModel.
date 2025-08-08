# This scripts runs a null model that preserves the number of plant species across null triads generated for the same site.
# The null model is focused on the plant-flower visitor interactions sampled on sites combining three different habitats (triads).
# In this null model variation, the number of interaction events to sample is preserved.

task_id <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(task_id[1])

# set number of repetitions
reps <- 1

# load libraries
lib2load <- c("permute","vegan", "statnet.common", "vctrs",
              "network", "sna", "bipartite", "dplyr", "tibble", "ape", "tictoc",
              "usedist", "geometry", "FD")
lapply(lib2load, function(x){library(x, character.only = TRUE, quietly = F)})
tic()
print('Packages loaded')

# define working directory
FXN_DIR <- "../functions"
DAT_DIR <- "../../data" # change path to your working directory if need be
OUT_DIR <- "../outputs"

# load functions
source(paste(FXN_DIR, "PD.R", sep = "/")) # function to calculate phylogenetic diversity
source(paste(FXN_DIR, "functional_diversity.R", sep = "/")) # functions for functional diversity
print('Functions sourced')

# load data
site_data <- read.table(paste(DAT_DIR, "site_level_data.txt", sep = "/"), header = T, sep = "\t", stringsAsFactors = FALSE)
allwebdata <- read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
pl_tree <- read.tree(paste(DAT_DIR, "plant_phylogeny/DaPhnE_01.tre", sep = "/")) # phylogenetic tree for plants
pl_solution <- read.table(paste(DAT_DIR, "plant_phylogeny/plant_solution.txt", sep = "/"), header = TRUE, sep = "\t", na.strings = "") # solutions for some plants that are not identified in the phylogenetic tree
print('Data loaded')

# create a flower-visitor metaweb for all sites
FV_metaweb <- allwebdata %>% filter(Web == "PL_FV")
FV_metaweb_mat <- frame2webs(FV_metaweb, varnames=c("Lower_Taxon", "Upper_Taxon", "Number"))[[1]] # dim: 150 524
dist_FV_metaweb <- bespokedist(FV_metaweb_mat, "sumto1", "bray") # the dissimilarity of diets within the upper guild
FV_metaweb <- FV_metaweb %>% select(-Habitat) %>% group_by(Site, Lower_Taxon, Upper_Taxon) %>% summarise(n_int = length(Upper_Taxon))

# create Flower-Visitor meta web for interactions within each habitat monad
FV_monads <- allwebdata %>% filter((Web == "PL_FV") & (M_D_T == "monad")) %>%
                            group_by(Habitat, Lower_Taxon, Upper_Taxon) %>%
                            summarise(n_int = length(Upper_Taxon))

# triad webs
FV_triads <- allwebdata %>% filter((Web == "PL_FV") & (M_D_T == "triad")) %>%
                            group_by(Site, Habitat, Lower_Taxon, Upper_Taxon) %>%
                            summarise(no_int = length(Upper_Taxon))

print('triads created')

# select triad site names in alphabetical order
triad_names <- sort(unique(FV_triads$Site))
N_triads <- length(triad_names)
N_hab_triads <- 3 # this is trivial, but avoids hard variables in loops

# create data frames to save network properties under investigation
nbdim_FD <- 10
FD_variables <- c("FRic_occ", "FEve_occ", "FDis_occ", "FDiv_occ", "FEve_abund", "FDis_abund", "FDiv_abund")
null_metrics <- data.frame(site = rep(triad_names, each = reps), rep = rep(1:reps, N_triads),
                           pl_div = NA, pl_phyl_div = NA, conn = NA, int_even = NA)
for (dim in nbdim_FD){
  null_metrics[, paste(FD_variables, dim, sep = "")] <- NA
}


# select the habitat list for each triad
for (Tr in 1:N_triads){
  site_web <- subset(FV_triads, Site == triad_names[Tr])
  T_habitats <- unique(site_web$Habitat)
  for (r in 1:reps){
    # for each habitat calculate the number of interaction events and select 
    # the same number from the meta-monad web
    # then aggregate the 3 habitat interaction lists into 1 null triad web
    # extract metric of interest and repeat
    null_triad <- tribble(~Habitat, ~Lower_Taxon, ~Upper_Taxon, ~n_int)
    for (H in 1:N_hab_triads){
      # print(c(Tr, r, H))
      # subset triad web to each habitat
      habitat_web <- subset(site_web, Habitat == T_habitats[H]) 
      
      # number of interaction events in 1 triad habitat web
      N_TH <- sum(habitat_web$no_int)
      
      # vector for each plant in 1 triad habitat web
      pl_species <- sort(unique(habitat_web$Lower_Taxon))
      
      # number of plant species
      N_TH_pl <- length(pl_species)
      
      # subset meta_monad to 1 habitat
      meta_monad_H <- subset(FV_monads, Habitat == T_habitats[H]); pl_meta_monad_H <- sort(unique(meta_monad_H$Lower_Taxon))
      monads_H_list <- sort(unique(allwebdata$Site[(allwebdata$Habitat == T_habitats[H]) & (allwebdata$M_D_T == "monad")])) # the names of the corresponding monads
      # get plant abundances in these monads according to plant surveys
      pl_abund_monads_H <- c()
      for (monad in monads_H_list){
        file_name <- paste(DAT_DIR, "/plant_abundances/", monad, "_PlAbundInt.csv", sep = "")
        pl_abund_monad <- read.csv(file_name, header = TRUE, stringsAsFactors = FALSE)
        pl_abund_monads_H <- rbind(pl_abund_monads_H, pl_abund_monad)
      }
      # keep the number of floral units FU because we focus on pollination interactions
      pl_abund_monads_H <- pl_abund_monads_H %>% filter(PlSpecies %in% pl_meta_monad_H) %>% group_by(PlSpecies) %>% summarise(FU = sum(FU))
      pl_abund_monads_H$FU[pl_abund_monads_H$FU == 0] <- 1
      pl_abund_monads_H <- pl_abund_monads_H[order(pl_abund_monads_H$PlSpecies), ]
      
      # vector for each plant species in 1 meta monad habitat web
      meta_monad_pl <- unique(meta_monad_H$Lower_Taxon)
      
      # sample the meta monad habitat web for a random selection of plant species the 
      # same length as the number of plants in the corresponding triad habitat
      pl_sp_null_triad <- sample(meta_monad_pl, size = N_TH_pl, prob = pl_abund_monads_H$FU, replace = FALSE)
      # plants are selected with probability proportional to their abundances
      
      # subset the meta monad habitat web to this list of plant species 
      meta_monad_H_plsubset <- subset(meta_monad_H, Lower_Taxon %in% pl_sp_null_triad)
      
      # sample the subsetted meta monad habitat web for the same number of interactions 
      # as in the triad habitat web weighted by the interactions in the monads
      int2sample <- sample(nrow(meta_monad_H_plsubset), size = N_TH, prob = meta_monad_H_plsubset$n_int, replace = TRUE)
      null_subtriad_H <- subset(meta_monad_H_plsubset[int2sample, ], select = c(Habitat, Lower_Taxon, Upper_Taxon))
      if (nrow(null_subtriad_H) != N_TH){
        stop("Problem with the number of sampled interaction events.")
      }
      
      # combine with other habitats
      null_triad <- bind_rows(null_triad, null_subtriad_H)
    }
    null_triad <- null_triad %>% select(-Habitat) %>% group_by(Lower_Taxon, Upper_Taxon) %>% summarise(n_int = length(Upper_Taxon))
    mat_null_triad <- xtabs(n_int ~ Lower_Taxon + Upper_Taxon, data = null_triad) # convert into a matrix
    if (sum(mat_null_triad) != sum(null_triad$n_int)){
      stop("Dataframe to matrix conversion went wrong.")
    }
    # calculate and save some network properties
    null_triad_metrics <- networklevel(mat_null_triad, index = c("connectance", "interaction evenness"), weighted = TRUE)
    null_metrics[(null_metrics$site == triad_names[Tr]) & (null_metrics$rep == r), c("conn", "int_even")] <- null_triad_metrics[c("connectance", "interaction evenness")]
    
    # diversity within the plant community
    pl_null_triad <- sort(unique(null_triad$Lower_Taxon))
    null_metrics$pl_div[(null_metrics$site == triad_names[Tr]) & (null_metrics$rep == r)] <- length(pl_null_triad)
    null_metrics$pl_phyl_div[(null_metrics$site == triad_names[Tr]) & (null_metrics$rep == r)] <- PD(pl_null_triad, pl_solution, pl_tree)[[1]]$mean_PD
    
    # functional diversity within the upper guild
    # create a data frame mixing all observed sites and the null triad - this ensures we can use the same trait space for all FD analyses
    null_triad <- bind_cols(data.frame(Site = rep("Null_triad", nrow(null_triad))), null_triad)
    all_int_null <- bind_rows(FV_metaweb, null_triad)
    null_triad_FD <- func_div_multsites(all_int_null, dist_FV_metaweb, nbdim_FD); null_triad_FD <- null_triad_FD[null_triad_FD$Site == "Null_triad", ]
    # save outputs of the functional diversity analysis
    for (n in nbdim_FD){
      null_metrics[(null_metrics$rep == r) & (null_metrics$site == triad_names[Tr]), paste(FD_variables, n, sep = "")] <- null_triad_FD[null_triad_FD$NoDim == n, FD_variables]
    }
  }
}
# export null model results
write.csv(null_metrics, paste(OUT_DIR, paste0("null_model_2_metrics_", task_id, ".csv"), sep = "/"), quote = FALSE, row.names = FALSE)

elapsed_time <- toc(quiet = TRUE); elapsed_time <- elapsed_time$toc - elapsed_time$tic
# exec_time <- data.frame(nm = "nm_2", ID = task_id, time = elapsed_time)
# write.table(exec_time, file = "exec_time.csv", sep = ",", append = file.exists("exec_time.csv"), quote = FALSE, col.names = !file.exists("exec_time.csv"), row.names = FALSE)
