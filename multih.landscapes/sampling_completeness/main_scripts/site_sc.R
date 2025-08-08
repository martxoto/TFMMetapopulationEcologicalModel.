# This script draws rarefaction curves/sampling coverage curves, and explores how we can subsample dataset to the same level of sampling completeness.
# All analyses rely on transposing the concepts of species richness estimation to interaction richness.

rm(list = ls(all = TRUE))

# load libraries
library(iNEXT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggrepel)
library(cowplot)

# define working directory
DAT_DIR <- "../../data"
OUT_DIR <- "../outputs"

# load data
all_obs <- read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_intSC <- unique(subset(all_obs, select = c(Site, M_D_T))); df_intSC$sc <- NA

for (site_name in df_intSC$Site){
  int_site <- all_obs %>% filter(Site == site_name) %>% group_by(Lower_Taxon, Upper_Taxon, Web) %>% count() %>% ungroup()
  out_inc <- iNEXT(int_site$n, datatype = "abundance", knots = 50, q = 0)
  df_intSC$sc[df_intSC$Site == site_name] <- out_inc$DataInfo$SC
}
df_intSC <- df_intSC[order(df_intSC$Site), ]

write.table(df_intSC, file = paste(OUT_DIR, "site_int-sampling-coverage.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)
