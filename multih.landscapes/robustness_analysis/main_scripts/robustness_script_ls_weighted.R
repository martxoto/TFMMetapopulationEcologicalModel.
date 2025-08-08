# This script performs a robustness analysis on the 30 sites prospected for the "landscape scale food webs" project.
# This robustness analysis includes rewiring, with interaction frequencies distributed in proportion to the abundance of the newly forage resource/prey.
# It also delete one species if one of its life stages cannot feed on anything.

rm(list=ls(all=TRUE))

library(abind)

# directories
FXN_DIR <- "../functions"
DAT_DIR <- "../../data"
OUT_DIR <- "../outputs"
  
# loading functions
fxn_list <- list.files(FXN_DIR)
for (fxn in fxn_list){
  source(paste(FXN_DIR, fxn, sep = "/"))
}

# loading the data
MDT_webs <- read.csv(paste(DAT_DIR, 'MDT_IntListAgg.csv', sep = "/"), header = TRUE, stringsAsFactors = FALSE)
MDT_PlAbund_files <- list.files(paste(DAT_DIR, "plant_abundances", sep = "/"), pattern = "PlAbundInt", recursive = TRUE)
MDT_PlAbund_files <- paste("plant_abundances", MDT_PlAbund_files, sep = "/")

# List of species from the different guilds
Guilds <- unique(c(MDT_webs$Lower_Guild, MDT_webs$Upper_Guild))
PL <- MDT_webs$Lower_Taxon[grep('PL_', MDT_webs$Web)]; PL <- unique(PL); PL <- sort(PL); PL <- as.character(PL)
LM <- c(MDT_webs$Lower_Taxon[grep('LM_', MDT_webs$Web)], MDT_webs$Upper_Taxon[grep('_LM', MDT_webs$Web)]); LM <- unique(LM); LM <- sort(LM); LM <- as.character(LM)
SF <- c(MDT_webs$Lower_Taxon[grep('SF_', MDT_webs$Web)], MDT_webs$Upper_Taxon[grep('_SF', MDT_webs$Web)]); SF <- unique(SF); SF <- sort(SF); SF <- as.character(SF)
FV <- MDT_webs$Upper_Taxon[grep('_FV', MDT_webs$Web)]; FV <- unique(FV); FV <- sort(FV); FV <- as.character(FV)
CP <- MDT_webs$Upper_Taxon[grep('_CP', MDT_webs$Web)]; CP <- unique(CP); CP <- sort(CP); CP <- as.character(CP)
P <- MDT_webs$Upper_Taxon[grep('_P', MDT_webs$Web)]; P <- unique(P); P <- sort(P); P <- as.character(P)
SFP <- MDT_webs$Upper_Taxon[grep('_SFP', MDT_webs$Web)]; SFP <- unique(SFP); SFP <- sort(SFP); SFP <- as.character(SFP)

sites <- unique(MDT_webs$Site); sites <- sort(sites) # list of sites
N_sites <- length(sites) # number of sites
web_types <- unique(subset(MDT_webs, select = c(Web, TL_lower_taxon, LS_lower_taxon, TL_upper_taxon, LS_upper_taxon))) # list of network types
web_types <- web_types[order(web_types$TL_lower_taxon),] # order interaction types
N_types <- nrow(web_types)

MDT_PlAbund <- vector("list", N_sites)
for (s in 1:N_sites){
  MDT_PlAbund[[s]] <- read.csv(paste(DAT_DIR, MDT_PlAbund_files[s], sep = "/"), header = TRUE, stringsAsFactors = FALSE)
}

# save plante relative abundances using only information about CP, and shaped as a matrix
# nrow(PL_BM[[i]] = no. of plant species), and ncol(PL_BM[[i]]) = number of habitats in site i
PL_BM <- vector("list", N_sites)
for (s in 1:N_sites){
  Hab_s <- sort(unique(MDT_PlAbund[[s]]$Habitat))
  NHab_s <- length(Hab_s)
  Pl_s <- sort(unique(MDT_webs$Lower_Taxon[(MDT_webs$Site == sites[s]) & (MDT_webs$Lower_Guild == "PL")]))
  NPl_s <- length(Pl_s)
  PL_BM[[s]] <- matrix(0, NPl_s, NHab_s, dimnames = list(Pl_s, Hab_s))
  TotCP_s <- vector("numeric", NHab_s); names(TotCP_s) <- Hab_s
  for (h in Hab_s){
    TotCP_s[h] <- sum(MDT_PlAbund[[s]]$CP[MDT_PlAbund[[s]]$Habitat == h])
  }
  for (p in Pl_s){
    CP_ps <- MDT_PlAbund[[s]]$CP[MDT_PlAbund[[s]]$PlSpecies == p]
    names(CP_ps) <- MDT_PlAbund[[s]]$Habitat[MDT_PlAbund[[s]]$PlSpecies == p]
    PL_BM[[s]][p, ] <- CP_ps/TotCP_s
  }
}

N_rep <- 500

setwd(OUT_DIR)

for (s1 in N_sites:1){
  WEB_s1 <- subset(MDT_webs, Site == sites[s1])
  
  # creates the plant list of site i, in alphabetical order
  PLANTS_s1 <- WEB_s1$Lower_Taxon[grep('PL_', WEB_s1$Web)]; PLANTS_s1 <- unique(PLANTS_s1); PLANTS_s1 <- sort(PLANTS_s1)
  N_PL_s1 <- length(PLANTS_s1)
  PL_BM_s1 <- PL_BM[[s1]]
  PL_C_s1 <- 1/ncol(PL_BM_s1)*rowSums(PL_BM_s1) # plant commonness based on local abundances

  # creates the insect list of site i, in alphabetical order
  INSECTS_s1 <- c(WEB_s1$Upper_Taxon, WEB_s1$Lower_Taxon[-grep('PL_', WEB_s1$Web)]); INSECTS_s1 <- unique(INSECTS_s1); INSECTS_s1 <- sort(INSECTS_s1)
  N_INS_s1 <- length(INSECTS_s1)
  
  # species list  & insects' life stages
  LS_s1 <- array(0, dim = c(N_INS_s1, nrow(web_types), 2), dimnames = list(INSECTS_s1, web_types$Web, c("larva", "adult")))
  SP_s1 <- matrix(0, nrow = length(c(PLANTS_s1, INSECTS_s1)), ncol = length(Guilds))
  rownames(SP_s1) <- c(PLANTS_s1, INSECTS_s1); colnames(SP_s1) <- Guilds
  for (i1 in INSECTS_s1){
    # to which guild insect i1 belongs to?
    LGtypes_i1 <- subset(WEB_s1, select = c(Web, Lower_Guild, LS_lower_taxon), WEB_s1$Lower_Taxon == i1); LGtypes_i1 <- unique(LGtypes_i1)
    UGtypes_i1 <- subset(WEB_s1, select = c(Web, Upper_Guild, LS_upper_taxon), WEB_s1$Upper_Taxon == i1); UGtypes_i1 <- unique(UGtypes_i1)
    types_i1 <- data.frame(Web = c(LGtypes_i1$Web, UGtypes_i1$Web), Guild = c(LGtypes_i1$Lower_Guild, UGtypes_i1$Upper_Guild),
                           Life_Stage = c(LGtypes_i1$LS_lower_taxon, UGtypes_i1$LS_upper_taxon))
    types_i1 <- unique(types_i1)
    SP_s1[i1, as.character(types_i1$Guild)] <- 1
    for (k in 1:nrow(types_i1)){
      LS_s1[i1, as.character(types_i1$Web[k]), as.character(types_i1$Life_Stage[k])] <- 1
    }
  }
  SP_s1[PLANTS_s1, 'PL'] <- 1
  
  # creates a N_PL_s1*N_INS_s1*nrow(web_types) array, each slice being the adjacency matrix between species (A_ijk = 1 if species j feeds on i with interaction type k) in a given interaction type
  A_s1 <- array(0, dim = c(N_PL_s1+N_INS_s1, N_PL_s1+N_INS_s1, N_types), dimnames = list(c(PLANTS_s1, INSECTS_s1), c(PLANTS_s1, INSECTS_s1), web_types$Web))
  N_INTi <- nrow(WEB_s1) # number of interactions in WEB_s1
  for (j1 in 1:N_INTi){
    LTj1 <- WEB_s1$Lower_Taxon[j1]; UTj1 <- WEB_s1$Upper_Taxon[j1]
    TYPEj1 <- WEB_s1$Web[j1]
    A_s1[LTj1, UTj1, TYPEj1] <- WEB_s1$Freq[j1]
  }

  # RobustSite_s1 <- data.frame(Site = character(), MDT = character(), SimNo = numeric(), ExtThrs = numeric(), FlexThrs = numeric(),
  #                             RMComm = numeric(), RLComm = numeric(), RRand = numeric())

  FileConn <- file(paste0(sites[s1], "_Robustness.csv"))
  writeLines(paste(c("Site", "MDT", "SimNo", "ExtThrs", "FlexThrs", "RMComm", "RLComm", "RRand"), collapse = ","), FileConn)
  close(FileConn)
  
  FileConn <- file(paste0(sites[s1], "_ExtSeq_Random.csv"))
  writeLines(paste(c("SimNo", "ExtThrs", "FlexThrs", paste("Prop", 0:N_PL_s1, sep = "")), collapse = ","), FileConn)
  close(FileConn)
  
  FileConn <- file(paste0(sites[s1], "_ExtSeq_MostComm.csv"))
  writeLines(paste(c("SimNo", "ExtThrs", "FlexThrs", paste("Prop", 0:N_PL_s1, sep = "")), collapse = ","), FileConn)
  close(FileConn)
  
  FileConn <- file(paste0(sites[s1], "_ExtSeq_LeastComm.csv"))
  writeLines(paste(c("SimNo", "ExtThrs", "FlexThrs", paste("Prop", 0:N_PL_s1, sep = "")), collapse = ","), FileConn)
  close(FileConn)
  
  for (ext_thrs in c(0.25, 0.5, 0.75)){
    for (flex_thrs in c(0, 0.25, 0.5, 1)){
      for (j1 in 1:N_rep){
        print(c(sites[s1], ext_thrs, flex_thrs, j1))
        RRand_j1 <- robustness_rw_ls_weighted(A_s1, SP_s1, PL_C_s1, LS_s1, ext_thrs, flex_thrs, "random")
        RMComm_j1 <- robustness_rw_ls_weighted(A_s1, SP_s1, PL_C_s1, LS_s1, ext_thrs, flex_thrs, "mostcommon")
        RLComm_j1 <- robustness_rw_ls_weighted(A_s1, SP_s1, PL_C_s1, LS_s1, ext_thrs, flex_thrs, "leastcommon")
        # Robust_j1 <- data.frame(Site = sites[s1], MDT = unique(WEB_s1$M_D_T), SimNo = j1, ExtThrs = ext_thrs, FlexThrs = flex_thrs,
        #                         RMComm = RMComm_j1$robustness["All species"],
        #                         RLComm = RLComm_j1$robustness["All species"],
        #                         RRand = RRand_j1$robustness["All species"])
        Robust_j1 <- c(sites[s1], unique(WEB_s1$M_D_T), j1, ext_thrs, flex_thrs,
                       RMComm_j1$robustness["All species"], RLComm = RLComm_j1$robustness["All species"],
                       RRand = RRand_j1$robustness["All species"])
        write(Robust_j1, file = paste(sites[s1], "_Robustness.csv", sep = ""), sep = ",", append = TRUE, ncolumns = length(Robust_j1))
        write(c(j1, ext_thrs, flex_thrs, as.vector(RRand_j1$remaining_sp[1, ])), file = paste(sites[s1], "_ExtSeq_Random.csv", sep = ""), sep = ",", append = TRUE, ncolumns = 4 + N_PL_s1)
        write(c(j1, ext_thrs, flex_thrs, as.vector(RMComm_j1$remaining_sp[1, ])), file = paste(sites[s1], "_ExtSeq_MostComm.csv", sep = ""), sep = ",", append = TRUE, ncolumns = 4 + N_PL_s1)
        write(c(j1, ext_thrs, flex_thrs, as.vector(RLComm_j1$remaining_sp[1, ])), file = paste(sites[s1], "_ExtSeq_LeastComm.csv", sep = ""), sep = ",", append = TRUE, ncolumns = 4 + N_PL_s1)
      } 
    }
  }
  # write.table(RobustSite_s1, file = paste(sites[s1], "_Robustness.csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE)
}
