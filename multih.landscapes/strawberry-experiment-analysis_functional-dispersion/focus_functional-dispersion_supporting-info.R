
#########################################################################
# -------------------- Functional Diversity Analysis -------------------- 

# Author: Kate Pereira Maia
# Date: 06/2021

# Analysis for the SI (similar to performed for the main text), but using 
# only species with high sampling.

# Part 1: Interaction matrix - defines the interaction matrix to be used (
# sensitivity analysis on sampling completeness is done here).
# Part 2: Calculates distances - defines the distance matrix to be used.
# Part 3: Evaluates functional space (code from Maire et al. 2015 - CHECKED).
# Part 4: Calculates functional diversity.
# Part 5: Plots and Statistical analyses.

# ---------------------- Loading data and packages ----------------------

library(bipartite); library(FD); library(vegan) 

DAT_DIR <- "../data" # change path to your working directory if need be

source("./functions/remove_abundance.R")
source("./functions/bespoke_distance.R")
source("./functions/mSD.R")

allwebdata <- read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = T, sep = ",", quote = "")
site_data <- read.table(paste(DAT_DIR, "site_level_data.txt", sep = "/"), header = T, sep = "\t")
sampcomp <-  read.table(paste(DAT_DIR, "insects_sampling-completeness.txt", sep = "/"), header = T, sep = "\t", quote = "")[,c(1,2,6)]
colnames(sampcomp) <- c("species", "n.indiv", "scomp") # data on sampling completenness

# Subset flower-Visitor interactions
FVwebdata <- subset(allwebdata, allwebdata$Web == "PL_FV") 

# Site names ordered alphabetically
site_data$Site <- gsub(" ", "_", site_data$Site); site_data$Site[21] <- "Pembrey"
site_names <- sort(site_data$Site)

# --------------------- Part 1 - Interaction Matrix ---------------------

# Sensitivity analysis for SAMPLING COMPLETENESS - use 0.5 OR 0.7 #
#well_sampled <- sampcomp[sampcomp$n.indiv >= 10 & sampcomp$scomp >= 0.5, "species"]
well_sampled <- sampcomp[sampcomp$n.indiv >= 10 & sampcomp$scomp >= 0.7, "species"]
FVwebdata <- FVwebdata[FVwebdata$Upper_Taxon %in% well_sampled,]

# Combining sites to create a pooled interaction matrix (METAWEB)
FV_metaweb <- frame2webs(FVwebdata, varnames = c("Lower_Taxon", "Upper_Taxon", "Number"), emptylist = TRUE)[[1]]
dim(FV_metaweb)
sum(rowSums(FV_metaweb) == 0); sum(colSums(FV_metaweb) == 0)

# -------------------- Part 2 - Calculates distances --------------------

# function to calculate distance measures
relbray <- bespoke_distance(FV_metaweb, "sumto1", "bray")
dist_metaweb <- relbray; dim(as.matrix(dist_metaweb))

# ------- Part 3 - Evaluates functional space (Maire et al. 2015) -------

meanSD <- mSD(dist_metaweb)

# Plotting the fit of different number of dimensions
plot(1, type = "n", axes = T, xlim = c(1, 45), ylim = c(min(meanSD), max(meanSD)),
     xlab = "Number of axes", ylab = "Quality of functional space (mSD)")
x <- c(2:(length(meanSD) + 1))
points(meanSD ~ x, pch = 16)
lines(meanSD ~ x, lty = 1)

# --------------- Part 4 - Calculates functional diversity --------------

FV <- levels(factor(FVwebdata$Upper_Taxon)) # names of flower visitors
abund <- matrix(0, length(site_names), length(FV))
rownames(abund) <- site_names
colnames(abund) <- FV

# FV species abundance (abund) in each site
for (i in 1:length(site_names)) {
  
  s <- site_names[i]
  ins <- table(FVwebdata[FVwebdata$Site == s, "Upper_Taxon"]) # insect abund in s
  abund[match(s, rownames(abund)), match(names(ins), colnames(abund))] <- ins
  
}

# Functional diversity metrics using 10 dimensions
ab_10D <- dbFD(dist_metaweb, abund, corr = "cailliez", stand.FRic = T, m = 10)

# --------------- Part 5 - Plots and Statistical analyses ---------------

all(site_names == site_data$Site) # Ordering sites CHECKED

data10D <- cbind(site_data[, 1:2], ab_10D$FEve, ab_10D$FDiv, ab_10D$FDis)
colnames(data10D)[-c(1:2)] <- c("abFEve", "abFDiv", "abFDis")
data10D$MDT <- factor(data10D$MDT); levels(data10D$MDT) <- c("M", "D", "T")

# ------- Results with 10D -------#

par(mfrow = c(1,3), mar = c(3,4,3,1), oma = c(2,2,3,1))
boxplot(abFEve ~ MDT, data = data10D, ylab = "Abundance")
boxplot(abFDiv ~ MDT, data = data10D)
boxplot(abFDis ~ MDT, data = data10D)
mtext("Functional Space with 10 Dimensions", 3, outer = TRUE)

mAbEve <- aov(abFEve ~ MDT,data = data10D); summary(mAbEve) # .
mAbDiv <- aov(abFDiv ~ MDT,data = data10D); summary(mAbDiv)
mAbDis <- aov(abFDis ~ MDT,data = data10D); summary(mAbDis) # **
TukeyHSD(mAbDis, conf.level = 0.95)

# removing diads to match strawberry analysis
data <- data10D[data10D$MDT %in% c("M", "T"),]
data$MDT <- as.numeric(factor(data$MDT))
summary(aov(abFDis ~ MDT,data = data))
t.test(data$MDT, data$abFDis)

dev.off()
boxplot(abFDis ~ MDT, data = data10D, xlab = "Number of Habitats", 
        ylab = "Functional Diversity", ylim = c(0.15, 0.7))
text(1:3, 0.65, c("a", "ab", "b"))
