# Food web properties are checked for each site (monads, dyads, triads), includinig connectance, interaction evenness,
# functional diversity, and plant phylogenetic richness.

# load libraries
lib2load <- c("bipartite", "dplyr", "tibble", "ape", "tictoc")
lapply(lib2load, function(x){library(x, character.only = TRUE, quietly = F)})

# define working directory
FXN_DIR <- "../functions"
DAT_DIR <- "../../data"
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

site_names <- sort(unique(allwebdata$Site))

# create a flower-visitor metaweb for all sites
FV_metaweb <- allwebdata %>% filter(Web == "PL_FV")
FV_metaweb_mat <- frame2webs(FV_metaweb, varnames = c("Lower_Taxon", "Upper_Taxon", "Number"))[[1]] # dim: 150 524
dist_FV_metaweb <- bespokedist(FV_metaweb_mat, "sumto1", "bray") # the dissimilarity of diets within the upper guild

# Create a list of species found in a triad that is not found in the corresponding
# monad habitat. Some species are quite common indiating incomplete sampling is a
# factor more so than rare species
triad_only_sp_list <- tribble(~Site, ~Habitat, ~Species) # create the dataframe
triad_names <- sort(unique(site_data$Site[site_data$MDT == "monad"])); N_triads <- length(triad_names)

for (Tr in 1:N_triads){
  site_web <- subset(FV_metaweb, FV_metaweb$Site == triad_names[Tr])
  triad_habitats <- unique(site_web$Habitat)
  for (H in 1:3){
    meta_monad_hab <- subset(FV_metaweb, (Habitat == triad_habitats[H]) & (M_D_T == "monad")) # subset dataframe on monads interactions to habitat H
    meta_monad_sp <- unique(meta_monad_hab$Lower_Taxon) # list of taxa in the lower guild in that monad
    triad_hab <- subset(site_web, Habitat == triad_habitats[H]) # subset triad to habitat H
    triad_pl_sp <- unique(triad_hab$Lower_Taxon) # list of taxa in the lower guild in that sub-triad
    triad_only_sp <- setdiff(triad_pl_sp, meta_monad_sp) # list of taxa in the lower guild that are found only in the sub-triad
    if (length(triad_only_sp)>0){
      for (sp in 1:length(triad_only_sp)){
        triad_only_sp_list <- triad_only_sp_list %>% add_row(Site = triad_hab$Site[1], Habitat = triad_hab$Habitat[1], Species = triad_only_sp[sp])
      }
    }
  }
}

write.csv(triad_only_sp_list, paste(OUT_DIR, "triad_only_species_list.csv", sep = "/"), quote = FALSE, row.names = FALSE)

nbdim_FD <- c(5, 10)
FD_variables <- c("FRic_occ", "FEve_occ", "FDis_occ", "FDiv_occ", "FEve_abund", "FDis_abund", "FDiv_abund")
obs_metrics <- data.frame(site = site_names, pl_div = NA, pl_phyl_div = NA, conn = NA, int_even = NA)
for (dim in nbdim_FD){
  obs_metrics[, paste(FD_variables, dim, sep = "")] <- NA
}

for (site in site_names){
  site_web <- FV_metaweb %>% filter(Site == site) %>% group_by(Lower_Taxon, Upper_Taxon) %>% summarise(no_int = length(Upper_Taxon))
  # network properties
  mat_obs_site <- xtabs(no_int ~ Lower_Taxon + Upper_Taxon, data = site_web)
  obs_site_metrics <- networklevel(mat_obs_site, index = c("connectance", "interaction evenness"), weighted = TRUE)
  obs_metrics[(obs_metrics$site == site), c("conn", "int_even")] <- obs_site_metrics[c("connectance", "interaction evenness")]
  # diversity within the plant community
  pl_obs <- sort(unique(site_web$Lower_Taxon))
  obs_metrics$pl_div[obs_metrics$site == site] <- length(pl_obs)
  obs_metrics$pl_phyl_div[obs_metrics$site == site] <- PD(pl_obs, pl_solution, pl_tree)[[1]]$mean_PD
}

# functional diversity within the upper guild
obs_FD <- func_div_multsites(FV_metaweb, dist_FV_metaweb, nbdim_FD)
# save outputs of the functional diversity analysis
for (n in nbdim_FD){
  obs_metrics[match(obs_metrics$site, unique(obs_FD$Site)), paste0(FD_variables, n)] <- obs_FD[obs_FD$NoDim == n, FD_variables]
}

# export results for the observed
write.csv(obs_metrics, paste(OUT_DIR, "obs_metrics.csv", sep = "/"), quote = FALSE, row.names = FALSE)