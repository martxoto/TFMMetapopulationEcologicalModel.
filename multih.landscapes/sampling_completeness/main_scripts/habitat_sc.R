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
DAT_DIR <- "../../data" # change path to your working directory if need be
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

# Is the number of interaction events a good proxy of the level of sampling completeness across sites? This is very unlikely.
# Let's take an illustration in the MDT dataset.
triad_name <- "Bystock"; monad_name <- "Wyeswood_Common"
# Both sites have a portion of their area covered with grassland (whether total or a third).
grassland_int <- int_FV_obs %>% filter(((Site == monad_name) | (Site == triad_name)) & (Habitat == "grassland")) %>%
                             group_by(Site, Lower_Taxon, Upper_Taxon) %>% count() %>%
                             ungroup() %>% spread(Site, n)
# a data frame for interaction names
grassland_int_names <- grassland_int %>% select(Lower_Taxon, Upper_Taxon)
grassland_int_names$Int_ID <- paste("Int", 1:nrow(grassland_int_names), sep = "_")
# subset the dataframe for interaction in grassland subset of these two sites
grassland_int <- subset(grassland_int, select = c(Bystock, Wyeswood_Common)); rownames(grassland_int) <- grassland_int_names$Int_ID
grassland_int[is.na(grassland_int)] <- 0; grassland_int <- as.data.frame(grassland_int)
out_inc <- iNEXT(grassland_int, datatype = "abundance", knots = 50, q = 0)

estimateD(grassland_int, datatype = "abundance", base = "size") # with the same number of individuals caught
est_cov_eg <- estimateD(grassland_int, datatype = "abundance", base = "coverage") # with the same level of sampling completeness
est_cov_eg$m[(est_cov_eg$site == monad_name) & (est_cov_eg$order == 0)]

pdf(paste(OUT_DIR, "eg_SC_grassland.pdf", sep = "/"), width = 7, height = 5)
SC_eg_plot <- ggiNEXT(out_inc, type = 2)

obs_eg_SC <- out_inc$DataInfo
sites_eg_col <- c("black", triads_col[triad_name]); names(sites_eg_col) <- c(monad_name, triad_name)
SC_eg_plot <- SC_eg_plot + scale_shape_manual(values = sites_pch[c(monad_name, triad_name)]) + 
  scale_colour_manual(values = sites_eg_col) +
  scale_fill_manual(values = sites_eg_col) +
  theme_bw() + 
  geom_label_repel(data = obs_eg_SC, aes(x = n, y = SC, label = site, colour = site), size = 2) +
  theme(legend.position = "none") +
  ylim(0, 1) +
  geom_vline(data = obs_eg_SC, aes(xintercept = n, colour = site), linetype = "dashed") +
  geom_vline(xintercept = est_cov_eg$m[(est_cov_eg$site == monad_name) & (est_cov_eg$order == 0)]) +
  geom_hline(data = obs_eg_SC, aes(yintercept = SC, colour = site), linetype = "dashed")
# theme(legend.position = "bottom", legend.box = "horizontal")
# hacking ggiNEXT to change the line width (c.f. https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html)
SC_eg_plot <- ggplot_build(SC_eg_plot)
SC_eg_plot$data[[2]]$size <- 0.5
SC_eg_plot$data[[1]]$size <- 3
SC_eg_plot <- ggplot_gtable(SC_eg_plot)
grid.draw(SC_eg_plot)
dev.off()

# Let's repeat the above for each habitat
monadtriad_comp <- data.frame(monad = character(), triad = character(), habitat = character(), SC_monad = numeric(), SC_triad = numeric(), n2sample = numeric()) # create and initiate a dataframe collecting the number of interaction events to sample in the monad
for (hab in hab_list){
  # hab <- hab_list[1]
  monads_hab <- unlist(unique(int_FV_obs %>% filter((Habitat == hab) & (M_D_T == "monad")) %>% select(Site))) # monads with this habitat
  triads_hab <- unlist(unique(int_FV_obs %>% filter((Habitat == hab) & (M_D_T == "triad")) %>% select(Site))) # triads with this habitat

  nb_monads_hab <- length(monads_hab); nb_triads_hab <- length(triads_hab)
  nb_sites_hab <- nb_monads_hab + nb_triads_hab
  
  # Both sites have a portion of their area covered with grassland (whether total or a third).
  hab_int <- int_FV_obs %>% filter((M_D_T != "dyad") & (Habitat == hab)) %>%
                         group_by(Site, Lower_Taxon, Upper_Taxon) %>% count() %>%
                         ungroup() %>% spread(Site, n)
  
  # what if we merge the monads?
  # hab_int <- int_FV_obs %>% filter((M_D_T != "dyad") & (Habitat == hab))
  # hab_int$Site[hab_int$M_D_T == "monad"] <- "Meta_monad"
  # hab_int <- hab_int %>% group_by(Site, Lower_Taxon, Upper_Taxon) %>% count() %>%
  #                        ungroup() %>% spread(Site, n)
  
  # a data frame for interaction names
  hab_int_names <- hab_int %>% select(Lower_Taxon, Upper_Taxon)
  hab_int_names$Int_ID <- paste("Int", 1:nrow(hab_int_names), sep = "_")
  # subset the dataframe for interaction in hab subset of these two sites
  hab_int <- subset(hab_int, select = -c(Lower_Taxon, Upper_Taxon)); rownames(hab_int) <- hab_int_names$Int_ID
  hab_int[is.na(hab_int)] <- 0; hab_int <- as.data.frame(hab_int)
  hab_int <- hab_int[, c(monads_hab, triads_hab)]
  
  out_inc <- iNEXT(hab_int, datatype = "abundance", knots = 50, q = 0)

  # define colours for each site
  sites_hab_col <- c(rep("black", nb_monads_hab), triads_col[triads_hab])
  names(sites_hab_col) <- colnames(hab_int)
  sites_hab_col <- sites_hab_col[order(names(sites_hab_col))]
  
  # define symbols shape for each site
  sites_hab_pch <- sites_pch[c(monads_hab, triads_hab)]
  names(sites_hab_pch) <- colnames(hab_int)
  
  SC_hab_plot <- ggiNEXT(out_inc, type = 2)

  obs_SC <- out_inc$DataInfo
  obs_SC$col <- sites_hab_col[obs_SC$site]
  SC_hab_plot <- SC_hab_plot + scale_shape_manual(values = sites_hab_pch) + 
                               scale_colour_manual(values = sites_hab_col) +
                               scale_fill_manual(values = sites_hab_col) +
                               theme_bw() + 
                               geom_label_repel(data = obs_SC, aes(x = n, y = SC, label = site, colour = site), size = 2) +
                               theme(legend.position = "none") +
                               ggtitle(hab_names[hab]) + ylim(0, 1)
                               # theme(legend.position = "bottom", legend.box = "horizontal")
  # hacking ggiNEXT to change the line width (c.f. https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html)
  SC_hab_plot <- ggplot_build(SC_hab_plot)
  SC_hab_plot$data[[2]]$size <- 0.5
  SC_hab_plot$data[[1]]$size <- 3
  SC_hab_plot <- ggplot_gtable(SC_hab_plot)
  # grid.draw(SC_hab_plot)
  assign(paste("SC", hab, "plot", sep = "_"), SC_hab_plot)
  
  # to assess the number of interaction events to sample in each monad for each {triad, habitat}-tuple, we need 
  # the sampling completeness within the monad to be equal or superior to those in the {triad, habitat}-tuple
  SC_monads <- out_inc$DataInfo$SC[match(monads_hab, out_inc$DataInfo$site)] # sampling completeness in the monads
  # SC_monads <- out_inc$DataInfo$SC[out_inc$DataInfo$site == "Meta_monad"] # sampling completeness in the meta-monad
  for (triad in triads_hab){
    SC_triad <- out_inc$DataInfo$SC[out_inc$DataInfo$site == triad] # sampling completeness in the {triad, habitat}-tuple
    n_triad <- out_inc$DataInfo$n[out_inc$DataInfo$site == triad]
    monads2sample <- which(SC_monads >= SC_triad)
    monadtriad_hab_comp <- data.frame(monad = monads_hab, triad = triad, habitat = hab, SC_monad = SC_monads, SC_triad = SC_triad, n2sample = NA)
    
    if (length(monads2sample) == 0){
      # n_monads <- out_inc$DataInfo$n[match(monads_hab, out_inc$DataInfo$site)]
      # print(c(hab, triad, n_triad, n_monads))
    }
    else {
      for (monad in monads_hab[monads2sample]){
        estimates_coverage <- estimateD(hab_int[, c(monad, triad)], datatype = "abundance", base = "coverage")
        monadtriad_hab_comp$n2sample[monadtriad_hab_comp$monad == monad] <- estimates_coverage$m[(estimates_coverage$site == monad) & (estimates_coverage$order == 0)]
      }
    }
    monadtriad_comp <- rbind(monadtriad_comp, monadtriad_hab_comp)
  }
}
write.table(monadtriad_comp, file = paste(OUT_DIR, "nb_int2sample.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)

pdf(paste(OUT_DIR, "SC_hab_sites.pdf", sep = "/"), width = 10, height = 7.5)
plot_grid(SC_grassland_plot, SC_heathland_plot,
          SC_salt_marsh_plot, SC_sand_dune_plot,
          SC_scrub_plot, SC_woodland_plot, labels = toupper(letters[1:6]))
dev.off()
