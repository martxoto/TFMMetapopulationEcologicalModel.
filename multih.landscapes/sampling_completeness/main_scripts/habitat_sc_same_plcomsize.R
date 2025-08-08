# This script draws rarefaction curves/sampling coverage curves, and explores how we can subsample dataset to the same level of sampling completeness.
# All analyses rely on transposing the concepts of species richness estimation to interaction richness.

# Here, we specifically check how focusing on a subset of the plant community within the monads change their sampling completeness.

rm(list = ls(all = TRUE))

# number of times a monad is subsampled with a given number of plant species
reps <- 100

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
int_FV_obs <- all_obs %>% filter(Web == "PL_FV")

hab_list <- sort(unique(int_FV_obs$Habitat))
hab_names <- c("Grassland", "Heathland", "Salt marsh", "Sand dune", "Scrub", "Woodland")
names(hab_names) <- hab_list
site_list <- sort(unique(int_FV_obs$Site)); nb_sites <- length(site_list)
triad_list <- sort(unique(int_FV_obs$Site[int_FV_obs$M_D_T == "triad"]))
monad_list <- sort(unique(int_FV_obs$Site[int_FV_obs$M_D_T == "monad"]))

triads_col <- c(brewer.pal(length(triad_list)-2, "Dark2"), brewer.pal(3, "Set2")[1:2])
names(triads_col) <- triad_list
sites_pch <- 1:(length(monad_list) + length(triad_list))
names(sites_pch) <- c(monad_list, triad_list)

# create and initiate a dataframe collecting the number of interaction events to sample in the monad
monadtriad_comp <- data.frame(monad = character(), triad = character(), habitat = character(),
                              n_pl_monad = numeric(), n_pl_triad = numeric(),
                              SC_triad = numeric(),
                              SC_monad_mean = numeric(), SC_monad_sd = numeric(),
                              n2sample_mean = numeric(), n2sample_sd = numeric())

# hereafter, we calculate the level of sampling completeness of each monads when subsampled with a given number of plant species
for (monad in monad_list){ # for each monad
  int_monad <- int_FV_obs %>% filter(Site == monad)%>% group_by(Site, Lower_Taxon, Upper_Taxon) %>% count()
  pl_monad <- sort(unique(int_monad$Lower_Taxon)); N_pl_monad <- length(pl_monad)
  hab_monad <- unique(int_FV_obs$Habitat[int_FV_obs$Site == monad])
  hab_triads_list <- unlist(unique(int_FV_obs %>% filter((Habitat == hab_monad) & (M_D_T == "triad")) %>% select(Site)))
  
  for (triad in hab_triads_list){
    
    int_triad <- int_FV_obs %>% filter((Site == triad) & (Habitat == hab_monad)) %>% group_by(Site, Lower_Taxon, Upper_Taxon) %>% count()
    pl_triad <- sort(unique(int_triad$Lower_Taxon)); N_pl_triad <- length(pl_triad)
    
    # create and initiate a line recording how the level of sampling completeness varies with the number of plants sampled
    # it also records how the number of interaction events to preserve varies
    monadtriad_hab_comp <- data.frame(monad = monad, triad = triad, habitat = hab_monad,
                                      n_pl_monad = N_pl_monad, n_pl_triad = N_pl_triad,
                                      SC_triad = NA,
                                      SC_monad_mean = NA, SC_monad_sd = NA,
                                      n2sample_mean = NA, n2sample_sd = NA)
    
    if (N_pl_triad <= N_pl_monad){
      SC_monad <- vector("numeric", reps); n2sample <- vector("numeric", reps)
      for (r in 1:reps){
        print(c(monad, triad, r))
        
        N_int_monad_sub <- 1
        # the following while loop ensures we sample networks with more than one link, otherwise estimating the number of interaction
        # events to sample to preserve sampling completeness is not possible
        while (N_int_monad_sub == 1){
          # vector for each plant species in 1 meta monad habitat web
          sub_pl_monad <- sample(pl_monad, N_pl_triad, replace = FALSE)
          
          # subset the meta monad habitat web to this list of plant species 
          int_monad_sub <- subset(int_monad, Lower_Taxon %in% sub_pl_monad)
          N_int_monad_sub <- nrow(int_monad_sub)
        }
        
        # bind the dataframes descripting the frequency of interaction {Lower_Taxon ; Upper_Taxon} in the monad and the triad (columns are named after sites)
        int_monadtriad <- int_monad_sub %>% bind_rows(int_triad) %>% spread(Site, n) 
        
        # a data frame for interaction names
        monadtriad_int_names <- int_monadtriad %>% select(Lower_Taxon, Upper_Taxon)
        monadtriad_int_names$Int_ID <- paste("Int", 1:nrow(monadtriad_int_names), sep = "_")
        
        # subset the dataframe for interaction in habitat hab subset of these two sites
        int_monadtriad <- subset(int_monadtriad, select = -c(Lower_Taxon, Upper_Taxon)); rownames(int_monadtriad) <- monadtriad_int_names$Int_ID
        int_monadtriad[is.na(int_monadtriad)] <- 0; int_monadtriad <- as.data.frame(int_monadtriad)
        out_inc <- iNEXT(int_monadtriad, datatype = "abundance", knots = 50, q = 0)
        SC_monad[r] <- out_inc$DataInfo$SC[out_inc$DataInfo$Assemblage == monad]
        if (r == 1){
          monadtriad_hab_comp$SC_triad <- out_inc$DataInfo$SC[out_inc$DataInfo$Assemblage == triad]
        }
        # estimated number of interactions events to sample in the subset of the monad
        est_cov <- estimateD(int_monadtriad, datatype = "abundance", base = "coverage") # with the same level of sampling completeness
        n2sample[r] <- est_cov$m[(est_cov$Assemblage == monad) & (est_cov$Order.q == 0)]
      }
      monadtriad_hab_comp$SC_monad_mean <- mean(SC_monad); monadtriad_hab_comp$SC_monad_sd <- sd(SC_monad)
      monadtriad_hab_comp$n2sample_mean <- mean(n2sample); monadtriad_hab_comp$n2sample_sd <- sd(n2sample)
    }
    monadtriad_comp <- rbind(monadtriad_comp, monadtriad_hab_comp)
  }
}

write.csv(monadtriad_comp, file = paste(OUT_DIR, "nb_int2sample_same_plcomsize.csv", sep = "/"), row.names = FALSE, quote = FALSE)
df_nb_int2sample <- read.csv(paste(OUT_DIR, "nb_int2sample.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE)

eg_monad <- "Wyeswood_Common"
eg_triads <- sort(as.character(monadtriad_comp$triad[monadtriad_comp$monad == eg_monad])); N_eg_monad <- length(eg_triads)

pdf(paste(OUT_DIR, "sc_n2sample_same_plcomsize_Wyeswood_Common.pdf", sep = "/"), width = 10, height = 5)
par(mfrow = c(1, 2))
# sampling completeness
plot(NA, NA, type = "n", xlim = c(0, N_eg_monad + 1), ylim = c(0, 1), xaxt = "n", xlab = NA, ylab = "Sampling completeness")
axis(1, at = 1:N_eg_monad, labels = FALSE)
text(x = 1:N_eg_monad, y = par("usr")[3] - 0.075, labels = eg_triads,
     xpd = NA, srt = 35, adj = 1)
triad_no <- 1
for (triad in eg_triads){
  index_monadtriad_tuple <- which((monadtriad_comp$monad == eg_monad) & (monadtriad_comp$triad == triad))
  segments(triad_no, monadtriad_comp$SC_monad_mean[index_monadtriad_tuple] - monadtriad_comp$SC_monad_sd[index_monadtriad_tuple],
           triad_no, monadtriad_comp$SC_monad_mean[index_monadtriad_tuple] + monadtriad_comp$SC_monad_sd[index_monadtriad_tuple],
           lwd = 10, col = adjustcolor("forestgreen", alpha.f = 0.5))
  points(triad_no, monadtriad_comp$SC_monad_mean[index_monadtriad_tuple], pch = 20)
  triad_no <- triad_no + 1
}
abline(h = unique(df_nb_int2sample$SC_monad[df_nb_int2sample$monad == eg_monad]), lty = 2)
legend("topleft", legend = "A", bty = "n")
# number of interaction events to sample
eg_monadtriad_comp <- subset(monadtriad_comp, monad == eg_monad)
plot(NA, NA, type = "n", xlim = c(0, N_eg_monad + 1), ylim = c(0, 700), xaxt = "n", xlab = NA, ylab = "No. of interaction events to sample")
axis(1, at = 1:N_eg_monad, labels = FALSE)
text(x = 1:N_eg_monad, y = par("usr")[3] - 50, labels = eg_triads,
     xpd = NA, srt = 35, adj = 1)
triad_no <- 1
for (triad in eg_triads){
  index_monadtriad_tuple <- which((monadtriad_comp$monad == eg_monad) & (monadtriad_comp$triad == triad))
  segments(triad_no, monadtriad_comp$n2sample_mean[index_monadtriad_tuple] - monadtriad_comp$n2sample_sd[index_monadtriad_tuple],
           triad_no, monadtriad_comp$n2sample_mean[index_monadtriad_tuple] + monadtriad_comp$n2sample_sd[index_monadtriad_tuple],
           lwd = 10, col = adjustcolor("forestgreen", alpha.f = 0.5))
  points(triad_no, monadtriad_comp$n2sample_mean[index_monadtriad_tuple], pch = 20)
  points(triad_no, df_nb_int2sample$n2sample[(df_nb_int2sample$monad == eg_monad) & (df_nb_int2sample$triad == triad)], pch = 17, col = "firebrick")
  triad_no <- triad_no + 1
}
legend("topleft", legend = "B", bty = "n")
dev.off()