
####################################################################################
#----------------------------- Bulding plant phylogeny -----------------------------
####################################################################################

#---------------------- Matching interaction data to phylogeny ---------------------

rm(list = ls())
DAT_DIR <- "../data"
OUT_DIR <- "outputs"

interac_data <- read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
interac_data <- interac_data[!interac_data$Web %in% c("CP_P", "LM_P"), ] # removes interactions with no plants
interac_data$Web <- factor(interac_data$Web)
interac_data$Lower_Taxon <- factor(interac_data$Lower_Taxon)
head(interac_data)

phyl_data <- read.table(paste(DAT_DIR, "plant_phylogeny/DaPhnE_02_Synonymy.txt", sep = "/"), header = TRUE, sep = "\t", quote = "")

# plant species per interaction type
plant_per_interac <- as.data.frame.matrix(table(interac_data$Lower_Taxon, interac_data$Web))
plant_list <- rownames(plant_per_interac)

# PL_CP
PL_CP_list <- rownames(plant_per_interac)[plant_per_interac[, "PL_CP"] > 0]
sum(PL_CP_list %in% phyl_data$Daphne_binomial) / length(PL_CP_list) # 53 of 76 in phylogeny
PL_CP_list[!PL_CP_list %in% phyl_data$Daphne_binomial]

# PL_FV
PL_FV_list <- rownames(plant_per_interac)[plant_per_interac[, "PL_FV"] > 0]
sum(PL_FV_list %in% phyl_data$Daphne_binomial) / length(PL_FV_list) # 132 of 150
PL_FV_list[!PL_FV_list %in% phyl_data$Daphne_binomial]

# PL_LM
PL_LM_list <- rownames(plant_per_interac)[plant_per_interac[, "PL_LM"] > 0]
sum(PL_LM_list %in% phyl_data$Daphne_binomial) / length(PL_LM_list) # 66 of 80
PL_LM_list[!PL_LM_list %in% phyl_data$Daphne_binomial]

# PL_SF
PL_SF_list <- rownames(plant_per_interac)[plant_per_interac[, "PL_SF"] > 0]
sum(PL_SF_list %in% phyl_data$Daphne_binomial) / length(PL_SF_list) # 17 of 20
PL_SF_list[!PL_SF_list %in% phyl_data$Daphne_binomial]

# PL_SFP
PL_SFP_list <- rownames(plant_per_interac)[plant_per_interac[, "PL_SFP"] > 0]
sum(PL_SFP_list %in% phyl_data$Daphne_binomial) / length(PL_SFP_list) # 16 of 18
PL_SFP_list[!PL_SFP_list %in% phyl_data$Daphne_binomial]

# Total
sum(plant_list %in% phyl_data$Daphne_binomial) / length(plant_list) # 169 of 212
miss_plant_list <- plant_list[!plant_list %in% phyl_data$Daphne_binomial] # MISS: sp in plant list not in phylogeny
ok_plant_list <- plant_list[plant_list %in% phyl_data$Daphne_binomial] # OK: sp in plant list & in phylogeny

data.frame(miss_plant_list)

write.table(data.frame(miss_plant_list), paste(OUT_DIR, "miss_plant_list.txt", sep = "/"), sep = "\t")

#---------------------- Fixing species names and finding solutions ---------------------

# Atention: Fixings/solutions for the miss_plant_list were found by eye in "plant_solution.xlsx"

miss_plant_list <- read.table(paste(DAT_DIR, "plant_phylogeny", "df_missing-plant-sp.txt", sep = "/"), header = TRUE, sep = "\t", na.strings = "")

interac_data$Phyl_Lower_Taxon <- interac_data$Lower_Taxon
all(levels(interac_data$Phyl_Lower_Taxon) == levels(interac_data$Lower_Taxon))

### By visual inspection, the 43 miss_plant can be divided into 4 groups:

### 1) Plant species with no possible solution (16)
miss_plant_list[miss_plant_list$solution == "N",]

### 2) Plant species with a simple solution (typos or synonyms) (13)
miss_plant_list[miss_plant_list$solution == "typo",]

### 3) Plant species which are the only representatives of the genus in the dataset (solution: sample representative) (5)
miss_plant_list[miss_plant_list$solution == "genus",]

### 4) Plant species which are not the only representatives of the genus in the dataset
# For these I searched for solutions in the phylogeny dataset (daphne) following the steps:
# a. Using daphne I created a dataframe with species of the problem genus which are not in 
# the interaction dataset - alternative daphne in plant_solution.xlsx
# b. Classified species into likely, not in SW and not in UK using GBIF
# c. Downloaded UK occurrences of likely species (GBIFDownl.R to download, results in 
# GBIFDownl_0037864-200221144449610.csv)
miss_plant_list[miss_plant_list$solution == "alternative",]

#--------------------------------------------------------------------------------------------------
# a. Using daphne to create dataframe with species of the problem genus which are not in 
# the interaction dataset - alternative daphne in plant_solution.xlsx
sp_to_average <- miss_plant_list[miss_plant_list$solution == "alternative",1]
sp_to_average <- unlist(lapply(strsplit(as.character(sp_to_average), "_"), function(x) x <- x[1]))

df <- phyl_data[phyl_data$Genus %in% sp_to_average & !phyl_data$Daphne_binomial %in% ok_plant_list,]
df <- df[,-c(3,5,6)]
# write.table(df, paste(OUT_DIR, "alternative_daphne.txt", sep = "/"), sep = "\t")
# Added to plant_solution.xlsx and used for tasks b and c.

#---------------------- Creating solutions for "alternative" with GBIF ---------------------

library(sp)
library(stringr)
library(dplyr)
library(parzer)

likely <- read.csv(paste(DAT_DIR, "plant_phylogeny", "likely_plant_species_search.csv", sep = "/")) # created from df
sites <- read.csv(paste(DAT_DIR, "sites_coords.csv", sep = "/"), header = TRUE, sep = ",")
gbif <- read.csv(paste(DAT_DIR, "plant_phylogeny", "GBIF_downl_likely-plant-species.csv", sep = "/"), header = TRUE, sep = "\t")
gbif <- gbif[,c(7:10,14,22:24,33,37)] # removing useless columns
dim(gbif); head(gbif)

sum(likely$Taxon.name %in% gbif$species) # 55 out of 61 likely species were downloaded
likely$Taxon.name[!likely$Taxon.name %in% unique(gbif$species)]

# Site coordinates: from DMS to DD
sites$fix_lat <- parse_lat(sites$Latitude)
sites$fix_lon <- parse_lon(sites$Longitude)

# "Extremes" os sites locations
S <- min(sites$fix_lat); N <- max(sites$fix_lat)
W <- min(sites$fix_lon); E <- max(sites$fix_lon)

# GBIF occurrences within sites "extremes"
gbif_lat <- which(gbif$decimalLatitude >= S & gbif$decimalLatitude <= N)
gbif_lon <- which(gbif$decimalLongitude >= W & gbif$decimalLongitude <= E)
gbif <- gbif[intersect(gbif_lat, gbif_lon),]
dim(gbif)

# Exploring match of likely and GBIF species 
sum(likely$Taxon.name %in% gbif$species) # now 48 out of 61 likely species remain in GBIF dataset
likely$Taxon.name[!likely$Taxon.name %in% gbif$species] # 13 species not within "extremes"
not_within <- likely$Taxon.name[!likely$Taxon.name %in% gbif$species]
not_within[!not_within %in% levels(gbif$verbatimScientificName)] # these 2 are trully not within
not_within <- not_within[!not_within %in% levels(gbif$verbatimScientificName)] 

not_requested <- levels(gbif$species)[!levels(gbif$species) %in% likely$Taxon.name]
table(factor(gbif[gbif$species %in% not_requested, "verbatimScientificName"]))
likely$Taxon.name[c(6,43)]

likely$Taxon.name <- as.character(likely$Taxon.name)
gbif$species <- factor(gbif$species)

# Creating dataframe of alternatives
altern_sp <- miss_plant_list[miss_plant_list$solution == "alternative",]
altern_sp$genus <- unlist(lapply(strsplit(as.character(altern_sp$miss_plant_list), "_"), function(x)x <- x[1]))

# Adding number of GBIF occurrences (species column) to likely
table_occ <- data.frame(table(gbif$species))
colnames(likely) <- colnames(table_occ)[1] <- "species"; table_occ$species <- as.character(table_occ$species)
likely <- left_join(likely, table_occ, by = "species")

# Adding number of GBIF occurrences (verbatim column) to likely
verbatim <- likely[is.na(likely$Freq), "species"]
table_occ <- data.frame(table(gbif[gbif$verbatimScientificName %in% verbatim, "verbatimScientificName"]))
colnames(table_occ)[1] <- "species"; table_occ$species <- as.character(table_occ$species)
likely <- left_join(likely, table_occ, by = "species")

likely[!is.na(likely$Freq.y) & likely$Freq.y != 0, "Freq.x"] <- likely[!is.na(likely$Freq.y) & likely$Freq.y != 0, "Freq.y"]
likely$genus <- unlist(lapply(strsplit(likely$species, " "), function(x)x <- x[1]))
likely <- likely[,c(1,4,2)]; colnames(likely)[3] <- "occ"
likely <- likely[!is.na(likely$occ),]

likely <- likely[with(likely, order(genus, -occ)),]
likely$factor <- NA

genus <- unique(likely$genus)
for (i in 1: length(genus)) {
  l <- likely[likely$genus == genus[i],]
  if (nrow(l) == 1) {
    likely[likely$genus == genus[i], "factor"] <- 1000
  } else {
    likely[likely$genus == genus[i], "factor"] <- round(l$occ/c(l$occ[-1],l$occ[1]), digits = 2)
  }
}

# To select the group of species we could use as alternatives I used the criteria:
# Ordered from most common to most rare, select up to the first to be twice as common ad the next
# (up to the first factor >= 2 for each genus)
# Coincidently, this also results in up to 3 alternative species per genus (except for Ranunculus)

cbind(no = 1:nrow(likely), likely)
# Ranunculus: it should be R. ficaria which is already in the dataset (syn. of Ficaria verna). So I chose next 3.
# Primula: it should be P. acaulis and vulgaris, but P. vulgaris is the binomial used for both inputs in daphne.
# So I only kept P. vulgaris.

likely <- likely[c(1:5,15:17,20,24,32,35:37),]

# No alternative for convolvulus in GBIF
sp_to_average[!sp_to_average %in% genus]
phyl_data[grep("Convolvulus", phyl_data$Input_binomial),]
interac_data[grep("Convolvulus" ,interac_data$Lower_Taxon),]
# From the 3 available Convolvulus in phyl_data, one is also in the interac_data (C. arvensis),
# one is not present in the UK (C. lineatus), only remaining C. cantabrica. 

likely <- rbind(data.frame(species = "Convolvulus cantabrica", genus = "Convolvulus", occ = 0, factor = 0), likely)

altern_sp$which <- as.character(altern_sp$which)
altern_print <- data.frame(miss_plant_list = character(0), solution = character(0), which = character(0), genus = character(0))
for (i in 1:nrow(altern_sp)) {
  n <- which(altern_sp$genus[i] == likely$genus)
  alt <- altern_sp[i,]
  if (length(n) > 1) {
    alt <- as.data.frame(lapply(alt, rep, length(n)))
    alt$which <- gsub(" ", "_", as.character(likely[n, "species"]))
  } else {
    alt[, "which"] <- gsub(" ", "_", as.character(likely[match(altern_sp$genus[i], likely$genus), "species"]))
  }
  altern_print <- rbind(altern_print, alt)
}

table(altern_print$genus) == table(likely$genus)

miss_plant_list <- miss_plant_list[-which(miss_plant_list$solution == "alternative"),]
miss_plant_list <- rbind(miss_plant_list, altern_print[, -4])

write.table(miss_plant_list, paste(DAT_DIR, "plant_phylogeny", "plant_solution.txt", sep = "/"), sep = "\t", row.names = FALSE)
