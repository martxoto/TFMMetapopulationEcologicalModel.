# This script draws boxplots of network properties in triads following 2 null models and their variants.
# The first null model tests whether triads are the sum of their habitats.
# The second one is similar, but also preserves the size of the plant community.
# The null models can be constrained either on the number of interaction events, or the level of sampling completeness.

rm(list = ls(all = TRUE))

# load libraries
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(TeachingDemos)
library(pryr)
library(lme4)

# define working directories
DAT_DIR <- "../../data"
OUT_DIR <- "../outputs"


# load data
all_obs <- read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
int_FV_obs <- all_obs %>% filter(Web == "PL_FV")

# load null model outcomes
nm1_metrics <- read.csv(paste(OUT_DIR, "null_model_1_metrics_20200608.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE)
nm2_metrics <- read.csv(paste(OUT_DIR, "null_model_2_metrics_20200608.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE)
nm1_sc_metrics <- read.csv(paste(OUT_DIR, "null_model_1_sc_metrics_20200608.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE)
nm2_sc_metrics <- read.csv(paste(OUT_DIR, "null_model_2_sc_metrics_20200608.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE)

nm_types <- c("nm1", "nm2", "nm1_sc", "nm2_sc")
for (nm in nm_types){
  df_name <- paste(nm, "metrics", sep = "_")
  nm_df <- get(df_name)
  nm_df$site <- as.factor(nm_df$site)
  assign(df_name, nm_df)
}

triad_names <- sort(unique(nm1_metrics$site)); N_triads <- length(triad_names)
site_col <- brewer.pal(10, "Set3"); names(site_col) <- triad_names

# load network properties for the observed networks
obs_metrics <- read.csv(paste(OUT_DIR, "obs_metrics.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE)
obs_metrics <- obs_metrics[match(triad_names, obs_metrics$site), ]

nm_list <- c("nm1", "nm1_sc", "nm2", "nm2_sc")

# wheck whether differences between observed and null are significant
for (triad in triad_names){
   for (nm in nm_list){
      df_name <- paste(nm, "metrics", sep = "_")
      nm_metrics <- get(df_name)
      nm_metrics_triad <- subset(nm_metrics, site == triad)

      # connectance
      obs_metrics[obs_metrics$site == triad, paste("pl_phyl_div", nm, sep = "_")] <- (obs_metrics$pl_phyl_div[obs_metrics$site == triad] - mean(nm_metrics_triad$pl_phyl_div)) / sd(nm_metrics_triad$pl_phyl_div)
      # connectance
      obs_metrics[obs_metrics$site == triad, paste("conn", nm, sep = "_")] <- (obs_metrics$conn[obs_metrics$site == triad] - mean(nm_metrics_triad$conn)) / sd(nm_metrics_triad$conn)
      # interaction evenness
      obs_metrics[obs_metrics$site == triad, paste("int_even", nm, sep = "_")] <- (obs_metrics$int_even[obs_metrics$site == triad] - mean(nm_metrics_triad$int_even)) / sd(nm_metrics_triad$int_even)
   }
}

# weighted connectance
pdf(paste(OUT_DIR, "weighted_connectance.pdf", sep = "/"), width = 10, height = 7)
par(mfrow = c(2, 2), mar = c(5.5, 4, 1, 1))
# null model #1
boxplot(conn ~ factor(site, exclude = NULL), data = nm1_metrics, ylim = c(0, 0.3), outline = FALSE, xlab = NA,
        ylab = "Weighted connectance", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$conn, pch = 21, bg = site_col)
legend("topleft", legend = "A", bty = "n")
# null model #1 + preserving sampling completeness
boxplot(conn ~ factor(site, exclude = NULL), data = nm1_sc_metrics, ylim = c(0, 0.3), outline = FALSE, xlab = NA,
        ylab = "Weighted connectance", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$conn, pch = 21, bg = site_col)
legend("topleft", legend = "B", bty = "n")
# null model #2
boxplot(conn ~ factor(site, exclude = NULL), data = nm2_metrics, ylim = c(0, 0.3), outline = FALSE, xlab = NA,
        ylab = "Weighted connectance", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$conn, pch = 21, bg = site_col)
legend("topleft", legend = "C", bty = "n")
# null model #2 + preserving sampling completeness
boxplot(conn ~ factor(site, exclude = NULL), data = nm2_sc_metrics, ylim = c(0, 0.3), outline = FALSE, xlab = NA,
        ylab = "Weighted connectance", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$conn, pch = 21, bg = site_col)
legend("topleft", legend = "D", bty = "n")
dev.off()


# mean plant phylogenetic diversity
pdf(paste(OUT_DIR, "pl_phyl_div.pdf", sep = "/"), width = 10, height = 7)
par(mfrow = c(2, 2), mar = c(5.5, 4, 1, 1))
# null model #1
boxplot(pl_phyl_div ~ factor(site, exclude = NULL), data = nm1_metrics, ylim = c(300, 3000), outline = FALSE, xlab = NA,
        ylab = "Mean plant phyl. div.", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$pl_phyl_div, pch = 21, bg = site_col)
legend("topleft", legend = "A", bty = "n")
# null model #1 + preserving sampling completeness
boxplot(pl_phyl_div ~ factor(site, exclude = NULL), data = nm1_sc_metrics, ylim = c(300, 3000), outline = FALSE, xlab = NA,
        ylab = "Mean plant phyl. div.", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$pl_phyl_div, pch = 21, bg = site_col)
legend("topleft", legend = "B", bty = "n")
# null model #2
boxplot(pl_phyl_div ~ factor(site, exclude = NULL), data = nm2_metrics, ylim = c(300, 3000), outline = FALSE, xlab = NA,
        ylab = "Mean plant phyl. div.", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$pl_phyl_div, pch = 21, bg = site_col)
legend("topleft", legend = "C", bty = "n")
# null model #2 + preserving sampling completeness
boxplot(pl_phyl_div ~ factor(site, exclude = NULL), data = nm2_sc_metrics, ylim = c(300, 3000), outline = FALSE, xlab = NA,
        ylab = "Mean plant phyl. div.", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$pl_phyl_div, pch = 21, bg = site_col)
legend("topleft", legend = "D", bty = "n")
dev.off()

# interaction evenness
pdf(paste(OUT_DIR, "interaction_evenness.pdf", sep = "/"), width = 10, height = 7)
par(mfrow = c(2, 2), mar = c(5.5, 4, 1, 1))
# null model #1
boxplot(int_even ~ factor(site, exclude = NULL), data = nm1_metrics, ylim = c(0.3, 0.7), outline = FALSE, xlab = NA,
        ylab = "Interaction evenness", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$int_even, pch = 21, bg = site_col)
legend("topleft", legend = "A", bty = "n")
# null model #1 + preserving sampling completeness
boxplot(int_even ~ factor(site, exclude = NULL), data = nm1_sc_metrics, ylim = c(0.3, 0.7), outline = FALSE, xlab = NA,
        ylab = "Interaction evenness", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$int_even, pch = 21, bg = site_col)
legend("topleft", legend = "B", bty = "n")
# null model #2
boxplot(int_even ~ factor(site, exclude = NULL), data = nm2_metrics, ylim = c(0.3, 0.7), outline = FALSE, xlab = NA,
        ylab = "Interaction evenness", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$int_even, pch = 21, bg = site_col)
legend("topleft", legend = "C", bty = "n")
# null model #2 + preserving sampling completeness
boxplot(int_even ~ factor(site, exclude = NULL), data = nm2_sc_metrics, ylim = c(0.3, 0.7), outline = FALSE, xlab = NA,
        ylab = "Interaction evenness", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$int_even, pch = 21, bg = site_col)
legend("topleft", legend = "D", bty = "n")
dev.off()

# functional richness
pdf(paste(OUT_DIR, "functional_richness.pdf", sep = "/"), width = 10, height = 7)
par(mfrow = c(2, 2), mar = c(5.5, 4, 1, 1))
# null model #1
boxplot(FRic_occ10 ~ factor(site, exclude = NULL), data = nm1_metrics, ylim = c(0, 0.2), outline = FALSE, xlab = NA,
        ylab = "Interaction evenness", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$FRic_occ10, pch = 21, bg = site_col)
legend("topleft", legend = "A", bty = "n")
# null model #1 + preserving sampling completeness
boxplot(FRic_occ10 ~ factor(site, exclude = NULL), data = nm1_sc_metrics, ylim = c(0, 0.2), outline = FALSE, xlab = NA,
        ylab = "Interaction evenness", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$FRic_occ10, pch = 21, bg = site_col)
legend("topleft", legend = "B", bty = "n")
# null model #2
boxplot(FRic_occ10 ~ factor(site, exclude = NULL), data = nm2_metrics, ylim = c(0, 0.2), outline = FALSE, xlab = NA,
        ylab = "Interaction evenness", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$FRic_occ10, pch = 21, bg = site_col)
legend("topleft", legend = "C", bty = "n")
# null model #2 + preserving sampling completeness
boxplot(FRic_occ10 ~ factor(site, exclude = NULL), data = nm2_sc_metrics, ylim = c(0, 0.2), outline = FALSE, xlab = NA,
        ylab = "Interaction evenness", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$FRic_occ10, pch = 21, bg = site_col)
legend("topleft", legend = "D", bty = "n")
dev.off()

# Interaction complementarity
pdf(paste(OUT_DIR, "functional_dispersion.pdf", sep = "/"), width = 10, height = 7)
par(mfrow = c(2, 2), mar = c(5.5, 4, 1, 1))
# null model #1
boxplot(FDis_abund10 ~ factor(site, exclude = NULL), data = nm1_metrics, ylim = c(0.4, 0.7), outline = FALSE, xlab = NA,
        ylab = "Interaction complementarity", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$FDis_abund10, pch = 21, bg = site_col)
legend("topleft", legend = "A", bty = "n")
# null model #1 + preserving sampling completeness
boxplot(FDis_abund10 ~ factor(site, exclude = NULL), data = nm1_sc_metrics, ylim = c(0.4, 0.7), outline = FALSE, xlab = NA,
        ylab = "Interaction complementarity", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$FDis_abund10, pch = 21, bg = site_col)
legend("topleft", legend = "B", bty = "n")
# null model #2
boxplot(FDis_abund10 ~ factor(site, exclude = NULL), data = nm2_metrics, ylim = c(0.4, 0.7), outline = FALSE, xlab = NA,
        ylab = "Interaction complementarity", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$FDis_abund10, pch = 21, bg = site_col)
legend("topleft", legend = "C", bty = "n")
# null model #2 + preserving sampling completeness
boxplot(FDis_abund10 ~ factor(site, exclude = NULL), data = nm2_sc_metrics, ylim = c(0.4, 0.7), outline = FALSE, xlab = NA,
        ylab = "Interaction complementarity", xaxt = "n", col = site_col)
axis(1, at = 1:N_triads, labels = FALSE, tck = -diff(par("usr")[3:4])/25)
text(x = 1:N_triads, y = par("usr")[3] - diff(par("usr")[3:4])/20, labels = triad_names,
     xpd = NA, srt = 35, adj = 1)
points(obs_metrics$FDis_abund10, pch = 21, bg = site_col)
legend("topleft", legend = "D", bty = "n")
dev.off()

# relationship with plant phylogenetic diversity
pdf(paste(OUT_DIR, "networkmetrics_pldiv_nm1.pdf", sep = "/"), width = 10, height = 14)
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
plot(nm1_metrics$pl_phyl_div, nm1_metrics$conn, pch = 20, col = adjustcolor(site_col[nm1_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Connectance", log = "y", xlim = c(300, 3000), ylim = c(0.03, 0.6))
points(obs_metrics$pl_phyl_div, obs_metrics$conn, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "A", bty = "n")
plot(nm1_sc_metrics$pl_phyl_div, nm1_sc_metrics$conn, pch = 20, col = adjustcolor(site_col[nm1_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Connectance", log = "y", xlim = c(300, 3000), ylim = c(0.03, 0.6))
points(obs_metrics$pl_phyl_div, obs_metrics$conn, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "B", bty = "n")

plot(nm1_metrics$pl_phyl_div, nm1_metrics$int_even, pch = 20, col = adjustcolor(site_col[nm1_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction evenness", log = "y", xlim = c(300, 3000), ylim = c(0.4, 0.8))
points(obs_metrics$pl_phyl_div, obs_metrics$int_even, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "C", bty = "n")
plot(nm1_sc_metrics$pl_phyl_div, nm1_sc_metrics$int_even, pch = 20, col = adjustcolor(site_col[nm1_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction evenness", log = "y", xlim = c(300, 3000), ylim = c(0.4, 0.8))
points(obs_metrics$pl_phyl_div, obs_metrics$int_even, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "D", bty = "n")

plot(nm1_metrics$pl_phyl_div, nm1_metrics$FRic_occ10, pch = 20, col = adjustcolor(site_col[nm1_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Functional richness", log = "y", xlim = c(300, 3000))
points(obs_metrics$pl_phyl_div, obs_metrics$FRic_occ10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "E", bty = "n")
plot(nm1_sc_metrics$pl_phyl_div, nm1_sc_metrics$FRic_occ10, pch = 20, col = adjustcolor(site_col[nm1_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Functional richness", log = "y", xlim = c(300, 3000))
points(obs_metrics$pl_phyl_div, obs_metrics$FRic_occ10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "F", bty = "n")

plot(nm1_metrics$pl_phyl_div, nm1_metrics$FDis_abund10, pch = 20, col = adjustcolor(site_col[nm1_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction complementarity", log = "y", xlim = c(300, 3000))
points(obs_metrics$pl_phyl_div, obs_metrics$FDis_abund10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "G", bty = "n")
plot(nm1_sc_metrics$pl_phyl_div, nm1_sc_metrics$FDis_abund10, pch = 20, col = adjustcolor(site_col[nm1_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction complementarity", log = "y", xlim = c(300, 3000))
points(obs_metrics$pl_phyl_div, obs_metrics$FDis_abund10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "H", bty = "n")
legend("bottomright", legend = triad_names, pch = 21, pt.bg = site_col, bty = "n", ncol = 2)

dev.off()

# controlling the number of plant species
pdf(paste(OUT_DIR, "networkmetrics_pldiv_nm2.pdf", sep = "/"), width = 10, height = 14)
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
plot(nm2_metrics$pl_phyl_div, nm2_metrics$conn, pch = 20, col = adjustcolor(site_col[nm2_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Connectance", log = "y", xlim = c(10, 2500))
points(obs_metrics$pl_phyl_div, obs_metrics$conn, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "A", bty = "n")
plot(nm2_sc_metrics$pl_phyl_div, nm2_sc_metrics$conn, pch = 20, col = adjustcolor(site_col[nm2_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Connectance", log = "y", xlim = c(10, 2500))
points(obs_metrics$pl_phyl_div, obs_metrics$conn, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "B", bty = "n")

plot(nm2_metrics$pl_phyl_div, nm2_metrics$int_even, pch = 20, col = adjustcolor(site_col[nm2_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction evenness", log = "y", xlim = c(10, 2500))
points(obs_metrics$pl_phyl_div, obs_metrics$int_even, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "C", bty = "n")
plot(nm2_sc_metrics$pl_phyl_div, nm2_sc_metrics$int_even, pch = 20, col = adjustcolor(site_col[nm2_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction evenness", log = "y", xlim = c(10, 2500))
points(obs_metrics$pl_phyl_div, obs_metrics$int_even, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "D", bty = "n")

plot(nm2_metrics$pl_phyl_div, nm2_metrics$FRic_occ10, pch = 20, col = adjustcolor(site_col[nm2_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Functional richness", log = "y", xlim = c(10, 2500))
points(obs_metrics$pl_phyl_div, obs_metrics$FRic_occ10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "E", bty = "n")
plot(nm2_sc_metrics$pl_phyl_div, nm2_sc_metrics$FRic_occ10, pch = 20, col = adjustcolor(site_col[nm2_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Functional richness", log = "y", xlim = c(10, 2500))
points(obs_metrics$pl_phyl_div, obs_metrics$FRic_occ10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "F", bty = "n")

plot(nm2_metrics$pl_phyl_div, nm2_metrics$FDis_abund10, pch = 20, col = adjustcolor(site_col[nm2_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction complementarity", log = "y", xlim = c(10, 2500))
points(obs_metrics$pl_phyl_div, obs_metrics$FDis_abund10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "G", bty = "n")
plot(nm2_sc_metrics$pl_phyl_div, nm2_sc_metrics$FDis_abund10, pch = 20, col = adjustcolor(site_col[nm2_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction complementarity", log = "y", xlim = c(10, 2500))
points(obs_metrics$pl_phyl_div, obs_metrics$FDis_abund10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "H", bty = "n")
legend("bottomright", legend = triad_names, pch = 21, pt.bg = site_col, bty = "n", ncol = 2)
dev.off()


# the following chunk of code creates the fifth figure displayed in the main text
# each panel combine a plot of a network metric for random triads against plant phylogenetic diversity
# with an inset of the distribution of these values for each site as individual boxplots + observed value
# Layout : 
# A | B
# C | D
# where A & B represent outputs for nm1, and C & D for nm2
# A & C are for interaction evenness and B & D for functional diversity (name to change?)

if (!dir.exists(as.character(Sys.Date()))){dir.create(as.character(Sys.Date()))}
setwd(as.character(Sys.Date()))

pdf(paste(OUT_DIR, "figure4_randtriads_plphyldiv.pdf", sep = "/"), width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
# for each null model
plot(nm1_metrics$pl_phyl_div, nm1_metrics$int_even, pch = 20, col = adjustcolor(site_col[nm1_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction evenness", xlim = c(500, 3000), ylim = c(0.2, 0.7), yaxt = "n")
axis(2, at = seq(0, 0.9, 0.1), labels = seq(0, 0.9, 0.1))
lm_nm <- lm(int_even ~ pl_phyl_div * site, data = nm1_metrics)
# adjr2 <- summary(lm_nm)$adj.r.squared
# legend("topright", legend = bquote(italic(R)^2 == .(format(adjr2, digits = 3))) , bty = "n")
b0 <- coefficients(lm_nm)["(Intercept)"] # intercept
s0 <- coefficients(lm_nm)["pl_phyl_div"] # slope
for (triad in triad_names){
  b <- b0 ; s <- s0
  # draw regression line for each site
  if (triad != triad_names[1]){
    if (summary(lm_nm)$coefficients[paste0("site", triad), 4] < 0.05){
      b <- b0 + coefficients(lm_nm)[paste0("site", triad)] # intercept
    }
    if (summary(lm_nm)$coefficients[paste0("pl_phyl_div:site", triad), 4] < 0.05){
      s <- s0 + coefficients(lm_nm)[paste0("pl_phyl_div:site", triad)] # slope 
    }
  }
  x <- range(nm1_metrics$pl_phyl_div[nm1_metrics$site == triad])
  segments(x[1], b + s * x[1], x[2], b + s * x[2], col = site_col[triad], lwd = 2)
}
lm_obs <- lm(int_even ~ pl_phyl_div, data = obs_metrics) ; summary(lm_obs)
points(obs_metrics$pl_phyl_div, obs_metrics$int_even, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "A", bty = "n")
triad_nice_names <- unlist(lapply(strsplit(triad_names, split = "_"), function(X){paste(X, collapse = " ")}))
legend("bottomleft", legend = triad_nice_names, pch = 20, col = site_col, ncol = 2, cex = 0.7)

a %<a-% {boxplot(int_even ~ factor(site, exclude = NULL), data = nm1_metrics, outline = FALSE, xlab = NA,
                 ylab = NA,  col = site_col, xaxt = "n", yaxt = "n", ylim = c(0.45, 0.7), xlim = c(0.45, 0.7),
                 at = obs_metrics$int_even, boxwex = 0.01, cex.lab = 0.7, horizontal = TRUE);
  abline(0, 1);
  axis(2, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  axis(1, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  mtext("Observed values", 2, line = 1, cex = 0.7) ; mtext("Null model values", 1, line = 1, cex = 0.7);
  legend("topleft", legend = toupper(letters[2]), bty = "n")
}
subplot(a, x = grconvertX(c(0.55, 0.98), from = 'npc'), y = grconvertY(c(0.15, 0.6), from = 'npc'))

plot(nm1_metrics$pl_phyl_div, nm1_metrics$FDis_abund10, pch = 20, col = adjustcolor(site_col[nm1_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction complementarity", xlim = c(500, 3000), ylim = c(0.35, 0.65))
lm_nm <- lm(FDis_abund10 ~ pl_phyl_div * site, data = nm1_metrics)
b0 <- coefficients(lm_nm)["(Intercept)"] # intercept
s0 <- coefficients(lm_nm)["pl_phyl_div"] # slope
for (triad in triad_names){
  b <- b0 ; s <- s0
  # draw regression line for each site
  if (triad != triad_names[1]){
    if (summary(lm_nm)$coefficients[paste0("site", triad), 4] < 0.05){
      b <- b0 + coefficients(lm_nm)[paste0("site", triad)] # intercept
    }
    if (summary(lm_nm)$coefficients[paste0("pl_phyl_div:site", triad), 4] < 0.05){
      s <- s0 + coefficients(lm_nm)[paste0("pl_phyl_div:site", triad)] # slope 
    }
  }
  x <- range(nm1_metrics$pl_phyl_div[nm1_metrics$site == triad])
  segments(x[1], b + s * x[1], x[2], b + s * x[2], col = site_col[triad], lwd = 2)
}
lm_obs <- lm(FDis_abund10 ~ pl_phyl_div, data = obs_metrics) ; summary(lm_obs)
points(obs_metrics$pl_phyl_div, obs_metrics$FDis_abund10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "C", bty = "n")

# alternative layout
a %<a-% {boxplot(FDis_abund10 ~ factor(site, exclude = NULL), data = nm1_metrics, outline = FALSE, xlab = NA,
                 ylab = NA,  col = site_col, xaxt = "n", yaxt = "n", ylim = c(0.5, 0.7), xlim = c(0.5, 0.7),
                 at = obs_metrics$int_even, boxwex = 0.01, cex.lab = 0.7, horizontal = TRUE);
  abline(0, 1);
  axis(2, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  axis(1, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  mtext("Observed values", 2, line = 1, cex = 0.7) ; mtext("Null model values", 1, line = 1, cex = 0.7);
  legend("topleft", legend = toupper(letters[4]), bty = "n")
}
subplot(a, x = grconvertX(c(0.55, 0.98), from = 'npc'), y = grconvertY(c(0.15, 0.6), from = 'npc'))

plot(nm2_metrics$pl_phyl_div, nm2_metrics$int_even, pch = 20, col = adjustcolor(site_col[nm2_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction evenness", xlim = c(350, 3000), ylim = c(0.2, 0.7))
lm_nm <- lm(int_even ~ pl_phyl_div * site, data = nm2_metrics)
b0 <- coefficients(lm_nm)["(Intercept)"] # intercept
s0 <- coefficients(lm_nm)["pl_phyl_div"] # slope
for (triad in triad_names){
  b <- b0 ; s <- s0
  # draw regression line for each site
  if (triad != triad_names[1]){
    if (summary(lm_nm)$coefficients[paste0("site", triad), 4] < 0.05){
      b <- b0 + coefficients(lm_nm)[paste0("site", triad)] # intercept
    }
    if (summary(lm_nm)$coefficients[paste0("pl_phyl_div:site", triad), 4] < 0.05){
      s <- s0 + coefficients(lm_nm)[paste0("pl_phyl_div:site", triad)] # slope 
    }
  }
  x <- range(nm2_metrics$pl_phyl_div[nm2_metrics$site == triad], na.rm = TRUE)
  segments(x[1], b + s * x[1], x[2], b + s * x[2], col = site_col[triad], lwd = 2)
}
lm_obs <- lm(int_even ~ pl_phyl_div, data = obs_metrics) ; summary(lm_obs)
points(obs_metrics$pl_phyl_div, obs_metrics$int_even, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "E", bty = "n")

a %<a-% {boxplot(int_even ~ factor(site, exclude = NULL), data = nm2_metrics, outline = FALSE, xlab = NA,
                 ylab = NA,  col = site_col, xaxt = "n", yaxt = "n", ylim = c(0.45, 0.7), xlim = c(0.45, 0.7),
                 at = obs_metrics$int_even, boxwex = 0.01, cex.lab = 0.7, horizontal = TRUE);
  abline(0, 1);
  axis(2, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  axis(1, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  mtext("Observed values", 2, line = 1, cex = 0.7) ; mtext("Null model values", 1, line = 1, cex = 0.7);
  legend("topleft", legend = toupper(letters[6]), bty = "n")
}
subplot(a, x = grconvertX(c(0.55, 0.98), from = 'npc'), y = grconvertY(c(0.15, 0.6), from = 'npc'))

plot(nm2_metrics$pl_phyl_div, nm2_metrics$FDis_abund10, pch = 20, col = adjustcolor(site_col[nm2_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction complementarity", xlim = c(350, 3000), ylim = c(0.35, 0.65))
lm_nm <- lm(FDis_abund10 ~ pl_phyl_div * site, data = nm2_metrics)
b0 <- coefficients(lm_nm)["(Intercept)"] # intercept
s0 <- coefficients(lm_nm)["pl_phyl_div"] # slope
for (triad in triad_names){
  b <- b0 ; s <- s0
  # draw regression line for each site
  if (triad != triad_names[1]){
    if (summary(lm_nm)$coefficients[paste0("site", triad), 4] < 0.05){
      b <- b0 + coefficients(lm_nm)[paste0("site", triad)] # intercept
    }
    if (summary(lm_nm)$coefficients[paste0("pl_phyl_div:site", triad), 4] < 0.05){
      s <- s0 + coefficients(lm_nm)[paste0("pl_phyl_div:site", triad)] # slope 
    }
  }
  x <- range(nm2_metrics$pl_phyl_div[nm2_metrics$site == triad], na.rm = TRUE)
  segments(x[1], b + s * x[1], x[2], b + s * x[2], col = site_col[triad], lwd = 2)
}
lm_obs <- lm(FDis_abund10 ~ pl_phyl_div, data = obs_metrics) ; summary(lm_obs)
points(obs_metrics$pl_phyl_div, obs_metrics$FDis_abund10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = "G", bty = "n")

a %<a-% {boxplot(FDis_abund10 ~ factor(site, exclude = NULL), data = nm2_metrics, outline = FALSE, xlab = NA,
                 ylab = NA,  col = site_col, xaxt = "n", yaxt = "n", ylim = c(0.5, 0.7), xlim = c(0.5, 0.7),
                 at = obs_metrics$int_even, boxwex = 0.01, cex.lab = 0.7, horizontal = TRUE);
  abline(0, 1);
  axis(2, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  axis(1, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  mtext("Observed values", 2, line = 1, cex = 0.7) ; mtext("Null model values", 1, line = 1, cex = 0.7);
  legend("topleft", legend = toupper(letters[8]), bty = "n")
}
subplot(a, x = grconvertX(c(0.55, 0.98), from = 'npc'), y = grconvertY(c(0.15, 0.6), from = 'npc'))
dev.off()

# the same but controlling for sampling completeness
pdf(paste(OUT_DIR, "fig18-ext-data_randtriads_sc_plphyldiv.pdf", sep = "/"), width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
# for each null model
plot(nm1_sc_metrics$pl_phyl_div, nm1_sc_metrics$int_even, pch = 20, col = adjustcolor(site_col[nm1_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction evenness", xlim = c(200, 2500), ylim = c(0.2, 0.7))
lm_nm <- lm(int_even ~ pl_phyl_div * site, data = nm1_sc_metrics)
adjr2 <- summary(lm_nm)$adj.r.squared
b0 <- coefficients(lm_nm)["(Intercept)"] # intercept
s0 <- coefficients(lm_nm)["pl_phyl_div"] # slope
triad_names_sc <- unique(nm1_sc_metrics$site[!is.na(nm1_sc_metrics$pl_phyl_div)])
for (triad in triad_names_sc){
  b <- b0 ; s <- s0
  # draw regression line for each site
  if (triad != triad_names_sc[1]){
    if (summary(lm_nm)$coefficients[paste0("site", triad), 4] < 0.05){
      b <- b0 + coefficients(lm_nm)[paste0("site", triad)] # intercept
    }
    if (summary(lm_nm)$coefficients[paste0("pl_phyl_div:site", triad), 4] < 0.05){
      s <- s0 + coefficients(lm_nm)[paste0("pl_phyl_div:site", triad)] # slope 
    }
  }
  x <- range(nm1_sc_metrics$pl_phyl_div[nm1_sc_metrics$site == triad], na.rm = TRUE)
  segments(x[1], b + s * x[1], x[2], b + s * x[2], col = site_col[triad], lwd = 2)
}
lm_obs <- lm(int_even ~ pl_phyl_div, data = obs_metrics) ; summary(lm_obs)
points(obs_metrics$pl_phyl_div, obs_metrics$int_even, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = toupper(letters[1]), bty = "n")

a %<a-% {boxplot(int_even ~ factor(site, exclude = NULL), data = nm1_sc_metrics, outline = FALSE, xlab = NA,
                 ylab = NA,  col = site_col, xaxt = "n", yaxt = "n", ylim = c(0.45, 0.7), xlim = c(0.45, 0.7),
                 at = obs_metrics$int_even, boxwex = 0.01, cex.lab = 0.7, horizontal = TRUE);
  abline(0, 1);
  axis(2, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  axis(1, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  mtext("Observed values", 2, line = 1, cex = 0.7) ; mtext("Null model values", 1, line = 1, cex = 0.7);
  legend("topleft", legend = toupper(letters[2]), bty = "n")
}
subplot(a, x = grconvertX(c(0.55, 0.98), from = 'npc'), y = grconvertY(c(0.15, 0.6), from = 'npc'))


plot(nm1_sc_metrics$pl_phyl_div, nm1_sc_metrics$FDis_abund10, pch = 20, col = adjustcolor(site_col[nm1_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction complementarity", xlim = c(200, 3000), ylim = c(0.25, 0.65))
lm_nm <- lm(FDis_abund10 ~ pl_phyl_div * site, data = nm1_sc_metrics)
adjr2 <- summary(lm_nm)$adj.r.squared
b0 <- coefficients(lm_nm)["(Intercept)"] # intercept
s0 <- coefficients(lm_nm)["pl_phyl_div"] # slope
for (triad in triad_names_sc){
  b <- b0 ; s <- s0
  # draw regression line for each site
  if (triad != triad_names_sc[1]){
    if (summary(lm_nm)$coefficients[paste0("site", triad), 4] < 0.05){
      b <- b0 + coefficients(lm_nm)[paste0("site", triad)] # intercept
    }
    if (summary(lm_nm)$coefficients[paste0("pl_phyl_div:site", triad), 4] < 0.05){
      s <- s0 + coefficients(lm_nm)[paste0("pl_phyl_div:site", triad)] # slope 
    }
  }
  x <- range(nm1_sc_metrics$pl_phyl_div[nm1_sc_metrics$site == triad], na.rm = TRUE)
  segments(x[1], b + s * x[1], x[2], b + s * x[2], col = site_col[triad], lwd = 2)
}
lm_obs <- lm(FDis_abund10 ~ pl_phyl_div, data = obs_metrics) ; summary(lm_obs)
points(obs_metrics$pl_phyl_div, obs_metrics$FDis_abund10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = toupper(letters[3]), bty = "n")

a %<a-% {boxplot(FDis_abund10 ~ factor(site, exclude = NULL), data = nm1_sc_metrics, outline = FALSE, xlab = NA,
                 ylab = NA,  col = site_col, xaxt = "n", yaxt = "n", ylim = c(0.45, 0.7), xlim = c(0.45, 0.7),
                 at = obs_metrics$int_even, boxwex = 0.01, cex.lab = 0.7, horizontal = TRUE);
  abline(0, 1);
  axis(2, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  axis(1, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  mtext("Observed values", 2, line = 1, cex = 0.7) ; mtext("Null model values", 1, line = 1, cex = 0.7);
  legend("topleft", legend = toupper(letters[4]), bty = "n")
}
subplot(a, x = grconvertX(c(0.55, 0.98), from = 'npc'), y = grconvertY(c(0.15, 0.6), from = 'npc'))

plot(nm2_sc_metrics$pl_phyl_div, nm2_sc_metrics$int_even, pch = 20, col = adjustcolor(site_col[nm2_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction evenness", xlim = c(10, 2200), ylim = c(-0.2, 0.95), yaxt = "n")
axis(2, at = seq(0, 0.9, 0.1), labels = seq(0, 0.9, 0.1))
lm_nm <- lm(int_even ~ pl_phyl_div * site, data = nm2_sc_metrics)
adjr2 <- summary(lm_nm)$adj.r.squared
b0 <- coefficients(lm_nm)["(Intercept)"] # intercept
s0 <- coefficients(lm_nm)["pl_phyl_div"] # slope
triad_names_sc <- unique(nm2_sc_metrics$site[!is.na(nm2_sc_metrics$FDis_abund10)])
for (triad in triad_names_sc){
  b <- b0 ; s <- s0
  # draw regression line for each site
  if (triad != triad_names_sc[1]){
    if (summary(lm_nm)$coefficients[paste0("site", triad), 4] < 0.05){
      b <- b0 + coefficients(lm_nm)[paste0("site", triad)] # intercept
    }
    if (summary(lm_nm)$coefficients[paste0("pl_phyl_div:site", triad), 4] < 0.05){
      s <- s0 + coefficients(lm_nm)[paste0("pl_phyl_div:site", triad)] # slope 
    }
  }
  x <- range(nm2_sc_metrics$pl_phyl_div[nm2_sc_metrics$site == triad], na.rm = TRUE)
  segments(x[1], b + s * x[1], x[2], b + s * x[2], col = site_col[triad], lwd = 2)
}
lm_obs <- lm(int_even ~ pl_phyl_div, data = obs_metrics) ; summary(lm_obs)
points(obs_metrics$pl_phyl_div, obs_metrics$int_even, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = toupper(letters[5]), bty = "n")
legend("topright", legend = triad_nice_names, pch = 20, col = site_col, ncol = 4, cex = 0.7)


a %<a-% {boxplot(int_even ~ factor(site, exclude = NULL), data = nm2_sc_metrics, outline = FALSE, xlab = NA,
                 ylab = NA,  col = site_col, xaxt = "n", yaxt = "n", ylim = c(0.2, 0.7), xlim = c(0.2, 0.7),
                 at = obs_metrics$int_even, boxwex = 0.01, cex.lab = 0.7, horizontal = TRUE);
  abline(0, 1);
  axis(2, at = seq(0.2, 0.7, by = 0.05), labels = seq(0.2, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  axis(1, at = seq(0.2, 0.7, by = 0.05), labels = seq(0.2, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  mtext("Observed values", 2, line = 1, cex = 0.7) ; mtext("Null model values", 1, line = 1, cex = 0.7);
  legend("topleft", legend = toupper(letters[6]), bty = "n")
}
subplot(a, x = grconvertX(c(0.55, 0.98), from = 'npc'), y = grconvertY(c(0.15, 0.6), from = 'npc'))

plot(nm2_sc_metrics$pl_phyl_div, nm2_sc_metrics$FDis_abund10, pch = 20, col = adjustcolor(site_col[nm2_sc_metrics$site], alpha.f = 0.2),
     xlab = "Plant phylogenetic diversity", ylab = "Interaction complementarity", xlim = c(10, 2200), ylim = c(0, 0.65))
lm_nm <- lm(FDis_abund10 ~ pl_phyl_div * site, data = nm2_sc_metrics)
adjr2 <- summary(lm_nm)$adj.r.squared
b0 <- coefficients(lm_nm)["(Intercept)"] # intercept
s0 <- coefficients(lm_nm)["pl_phyl_div"] # slope
triad_names_sc <- unique(nm2_sc_metrics$site[!is.na(nm2_sc_metrics$FDis_abund10)])
for (triad in triad_names_sc){
  b <- b0 ; s <- s0
  # draw regression line for each site
  if (triad != triad_names_sc[1]){
    if (summary(lm_nm)$coefficients[paste0("site", triad), 4] < 0.05){
      b <- b0 + coefficients(lm_nm)[paste0("site", triad)] # intercept
    }
    if (summary(lm_nm)$coefficients[paste0("pl_phyl_div:site", triad), 4] < 0.05){
      s <- s0 + coefficients(lm_nm)[paste0("pl_phyl_div:site", triad)] # slope 
    }
  }
  x <- range(nm2_sc_metrics$pl_phyl_div[nm2_sc_metrics$site == triad], na.rm = TRUE)
  segments(x[1], b + s * x[1], x[2], b + s * x[2], col = site_col[triad], lwd = 2)
}
lm_obs <- lm(FDis_abund10 ~ pl_phyl_div, data = obs_metrics) ; summary(lm_obs)
points(obs_metrics$pl_phyl_div, obs_metrics$FDis_abund10, pch = 21, bg = site_col[obs_metrics$site])
legend("topleft", legend = toupper(letters[7]), bty = "n")

a %<a-% {boxplot(FDis_abund10 ~ factor(site, exclude = NULL), data = nm2_sc_metrics, outline = FALSE, xlab = NA,
                 ylab = NA,  col = site_col, xaxt = "n", yaxt = "n", ylim = c(0.4, 0.7), xlim = c(0.4, 0.7),
                 at = obs_metrics$int_even, boxwex = 0.01, cex.lab = 0.7, horizontal = TRUE);
  abline(0, 1);
  axis(2, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  axis(1, at = seq(0.4, 0.7, by = 0.05), labels = seq(0.4, 0.7, by = 0.05), tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.3, 0));
  mtext("Observed values", 2, line = 1, cex = 0.7) ; mtext("Null model values", 1, line = 1, cex = 0.7);
  legend("topleft", legend = toupper(letters[8]), bty = "n")
}
subplot(a, x = grconvertX(c(0.55, 0.98), from = 'npc'), y = grconvertY(c(0.15, 0.6), from = 'npc'))

dev.off()

