
# +++++++++++++++++++++++++++++++++++++++++++
#### Tidied version for Nature, Mar 2024 ####

# Analysis for Suppl. Mat. Section 1:Variation associated with shared habitats

# Load required packages
library(vegan)
library(reshape2)

DAT_DIR <- "../data"

# 1. Variance partitioning to estimate how much variation in the plant 
#    communities among triads could be accounted for by study site versus 
#    habitat type.

# Load data and select data from 'triads'.
d1 <- read.csv(paste(DAT_DIR, "plant_abundances", "floral-units_per-site-hab.csv", sep = "/"))
vp.d <- d1[d1$Type == "tri", ]


# Results including all 30 triad sites (3 habitats x 10 locations)
vp1 <- varpart(vegdist(log1p(vp.d[, -c(1:4)]), method = "bray"), 
               data.frame(vp.d$Site), data.frame(vp.d$Habitat))

# Output table. Following the help page for varpart:
#    [a] = unique fraction explained by Site
#    [c] = unique fraction explained by Habitat
vp1


# Same analysis, but removing one row (habitat) that had zero floral 
#   abundance - hence generated a warning message. The results are very 
#   similar so the analysis including all 30 observations (vp1) so the 
#   full analysis is reported. 
vp.d2 <- vp.d[-3, ]
rowSums(vp.d2[, -c(1:4)])

vp2 <- varpart(vegdist(log1p(vp.d2[, -c(1:4)]), method = "bray"), 
               data.frame(vp.d2$Site), data.frame(vp.d2$Habitat))
vp2



# 2. Test for differences in beta-diversity between monads, dyads and triads.

# Load data and cast into 'wide' format
d1 <- read.csv(paste(DAT_DIR, "plant_abundances", "floral-units_per-site.csv", sep = "/"))
d2 <- dcast(d1, Site + Type ~ Plant_Species, fill = 0)


# Calculate mean Bray-Curtis dissimilarities between plant communities
#   within the three treatment types (monad, dyad, triad)
d.mon <- vegdist(log1p(d2[d2$Type == "mon", -c(1:2)]), method = "bray")
d.di <- vegdist(log1p(d2[d2$Type == "di", -c(1:2)]), method = "bray")
d.tri <- vegdist(log1p(d2[d2$Type == "tri", -c(1:2)]), method = "bray")

mean(d.mon)
mean(d.di)
mean(d.tri)


# Test for differences in dispersion (beta-diversity) between the three 
#   treatment types.
anova(betadisper(vegdist(log1p(d2[, -c(1:2)]), method = "bray"), 
                 group = d2$Type))



# +++++++++++++++++++++++++++++++++++++++++++

