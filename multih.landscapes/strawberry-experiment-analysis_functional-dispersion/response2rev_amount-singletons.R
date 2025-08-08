
# Author: Kate Maia
# Date: 10/07/23

# Analysis for reponse to reviewers, for 2nd ressubmission.

# ---------- Loading library and data ---------- #

library(tidyverse)
library(cowplot)

DAT_DIR <- "../data" # change path to your working directory if need be

allint <-   read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = T, sep = ",", quote = "")
sitel <- read.table(paste(DAT_DIR, "site_level_data.txt", sep = "/"), header = T, sep = "\t")

head(allint); table(allint$Web) # all interaction types

# ---------- Defining SINGLETON function ---------- #
# Counts the number and proportion of singletons (unique occurrences) in a given site

singlcount <- function(x) { 
  tb <- table(x$Upper_Taxon) == 1
  data.frame(n = sum(tb), prp = sum(tb)/length(tb))
}

# ---------- Calculating and plotting SINGLETON results ---------- #

# ----- ONLY FV WEBS

allintFV <- allint[allint$Web == "PL_FV",]
allintFV$Site <- factor(allintFV$Site)

allintFVl <- split(allintFV, ~Site) # splits df per site
all(unlist(lapply(allintFVl, function(x) length(unique(x$Site)))) == 1) # checks 1 site/df

singl_siteFV <- do.call(rbind, lapply(allintFVl, singlcount)) # singleton per site
colnames(singl_siteFV) <- paste("FV", colnames(singl_siteFV), sep = "_")
all(sitel$Site == gsub("_", " ", rownames(singl_siteFV))) # all same but Pembrey/Pembrey Beach
sitel <- cbind(sitel, singl_siteFV)

pa <- sitel %>% mutate(MDT = factor(MDT, levels = c(1, 2, 3))) %>% # Absolute number
  mutate(MDT = recode_factor(MDT, '1' = "Monad", '2' = "Dyad", '3' = "Tryad")) %>% 
  ggplot(aes(x = MDT, y = FV_n)) + geom_boxplot() +
  xlab("Treatment") + ylab("Number of singletlons (FV webs)") + theme_cowplot(); pa

pb <- sitel %>% mutate(MDT = factor(MDT, levels = c(1, 2, 3))) %>% # Proportion
  mutate(MDT = recode_factor(MDT, '1' = "Monad", '2' = "Dyad", '3' = "Tryad")) %>% 
  ggplot(aes(x = MDT, y = FV_prp)) + geom_boxplot() +
  xlab("Treatment") + ylab("Prop singletlons (FV webs)") + theme_cowplot(); pb

plot_grid(pa, pb, labels = c("a)", "b)"))

# ----- ALL WEBS

table(allint$Web) # all interaction types
allint$Site <- factor(allint$Site)

allintl <- split(allint, ~Site) # splits df per site
all(unlist(lapply(allintl, function(x) length(unique(x$Site)))) == 1)

singl_site <- do.call(rbind, lapply(allintl, singlcount)) # singleton per site
colnames(singl_site) <- paste("all", colnames(singl_site), sep = "_")
all(sitel$Site == gsub("_", " ", rownames(singl_site))) # all same, but Pembrey/Pembrey Beach
sitel <- cbind(sitel, singl_site)

pa <- sitel %>% mutate(MDT = factor(MDT, levels = c(1, 2, 3))) %>% 
  mutate(MDT = recode_factor(MDT, '1' = "Monad", '2' = "Dyad", '3' = "Tryad")) %>% 
  ggplot(aes(x = MDT, y = all_n)) + geom_boxplot() +
  xlab("Treatment") + ylab("Number of singletlons (all webs)") + theme_cowplot(); pa
  
pb <- sitel %>% mutate(MDT = factor(MDT, levels = c(1, 2, 3))) %>% 
  mutate(MDT = recode_factor(MDT, '1' = "Monad", '2' = "Dyad", '3' = "Tryad")) %>% 
  ggplot(aes(x = MDT, y = all_prp)) + geom_boxplot() +
  xlab("Treatment") + ylab("Prop singletlons (all webs)") + theme_cowplot(); pb

plot_grid(pa, pb, labels = c("a)", "b)"))
