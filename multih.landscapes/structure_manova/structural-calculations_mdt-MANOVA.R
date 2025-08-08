#libraries
library(lme4)
library(LMERConvenienceFunctions)
library(tidyverse)
library(ggplot2)
library(lattice)
library(grid)
library(gridBase)
library(plotrix)

DAT_DIR <- "../data" # change path to your working directory if need be
OUT_DIR <- "outputs"

#import and attach interaction data
allwebdata <-   read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = T, sep = ",", quote = "")
attach(allwebdata)

#site names
sites <- unique(sort(Site))

#MDT
MDT <- unique(M_D_T)

#calculate abundance of insects per site and save as txt file
allwebdata_no0 <-
  subset(allwebdata, allwebdata$Upper_Taxon_abundance != "")
allinsectspersite <-
  aggregate(allwebdata_no0$Upper_Taxon_abundance,
            by = list(allwebdata_no0$Site),
            length)
write.table(allinsectspersite, file = paste(OUT_DIR, "insects_total-number_per-site.txt", sep = "/"), sep ="\t")


#calculate species richness for all insects per site and save as txt file
number_of_insect_sp_persite <-
  aggregate(
    allwebdata_no0$Upper_Taxon_abundance,
    by = list(allwebdata_no0$Site),
    FUN = function(x)
      length(unique(x))
  )

write.table(number_of_insect_sp_persite, paste(OUT_DIR, "insects_total-sp-richness_per-site.txt", sep = "/"), sep = "\t")

detach(allwebdata)

site_data <- read.table(paste(DAT_DIR, "site_level_data.txt", sep = "/"), header = T, sep = "\t")
insects_persite <-
  as.data.frame(cbind(
    site_data$MDT,
    allinsectspersite$x,
    number_of_insect_sp_persite$x
  ))
colnames(insects_persite) <-
  c("MDT", "total_insect_sp_abund", "total_insect_sp_rich")

#subset MDT
m <- subset(insects_persite, insects_persite$MDT == 1)
d <- subset(insects_persite, insects_persite$MDT == 2)
t <- subset(insects_persite, insects_persite$MDT == 3)

#read table for floral diversity
floral_diversity <- read.csv(paste(DAT_DIR, "MDT_floral-diversity.csv", sep = "/"), header = T, sep = ",", quote = "")
attach(floral_diversity)

#plant names
plant_names <- unique(sort(floral_diversity$Plant_Species))

#no_of plants
no_plants <- length(plant_names)

#floral diversity by site
group_fd <- group_by(floral_diversity, Site, Plant_Species)
floral_units <- summarise(group_fd, sum(FU_Open, na.rm = T))
colnames(floral_units)[3] <- "FU_Open"
write.table(floral_units, paste(OUT_DIR, "floral-units_per-sites.txt", sep = "/"), sep = "\t")

group_fr <- group_by(floral_units, Site)
floral_richness_persite <- summarise(group_fr, length(Plant_Species))
floral_abundance_persite <- summarise(group_fr, sum(FU_Open))
detach(floral_diversity)


floral_rich_abund <-
  as.data.frame(
    cbind(
      site_data$MDT,
      floral_richness_persite$`length(Plant_Species)`,
      floral_abundance_persite$`sum(FU_Open)`
    )
  )
colnames(floral_rich_abund) <-
  c("MDT", "plant_rich", "openFU_abundance")

#subset MDT
p_m <- subset(floral_rich_abund, floral_rich_abund$MDT == 1)
p_d <- subset(floral_rich_abund, floral_rich_abund$MDT == 2)
p_t <- subset(floral_rich_abund, floral_rich_abund$MDT == 3)

###SPECIES EVENNESS CALCULATIONS

#libraries
library(vegan)

#read in data
allwebdata <-   read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = T, sep = ",", quote = "")
allwebdata_no0 <- subset(allwebdata, allwebdata$Upper_Taxon_abundance != "")
siteleveldata <- read.table(paste(DAT_DIR, "site_level_data.txt", sep = "/"), header = T, sep = "\t")

attach(allwebdata_no0)

#site names ordered alphabetically
site_names <- sort(unique(Site))
no_sites <- length(site_names)
MDT <- as.numeric(siteleveldata$MDT)

#species abundance table for all insect species at each site
sitespeciesmatrix <- vector("list", 30)
for (sn in 1:no_sites)
{
  sitedata <- subset(allwebdata, allwebdata_no0$Site == site_names[sn])
  sitedataspeciestable <-
    as.matrix(table(sitedata$Upper_Taxon_abundance))
  sitedataspecies <-
    subset(sitedataspeciestable, sitedataspeciestable[, 1] != 0)
  sitespeciesmatrix[[sn]] <- sitedataspecies
}

#species abundance table for all plant species at each site
p_sitespeciesmatrix <- vector("list", 30)
for (pn in 1:no_sites)
{
  p_sitedata <- subset(allwebdata, allwebdata_no0$Site == site_names[pn])
  p_only_sitedata <-
    subset(p_sitedata,
           p_sitedata$Web == c("PL_FV", "PL_LM", "PL_CP", "PL_SF"))
  p_sitedataspeciestable <-
    as.matrix(table(p_only_sitedata$Lower_Taxon))
  p_sitedataspecies <-
    subset(p_sitedataspeciestable, p_sitedataspeciestable[, 1] != 0)
  p_sitespeciesmatrix[[pn]] <- p_sitedataspecies
}

#calculate species evenness - Pilelou's #vegan package: diversity function to
#get shannon's diversity and then manually #calculate evenness as "Pielou's
#evenness J = H0/ log(S) is easily found as:> #J <- H/log(specnumber(BCI))"
J_insects <- matrix(nrow = 30, ncol = 3)
J_insects[, 1] <- as.numeric(site_names)
J_insects[, 2] <- as.numeric(MDT)
for (spe in 1:no_sites)
{
  H <- diversity(sitespeciesmatrix[[spe]])
  J <- H / log(length(sitespeciesmatrix[[spe]]))
  J_insects[spe, 3] <- as.numeric(J)
}

J_plants <- matrix(nrow = 30, ncol = 3)
J_plants[, 1] <- as.numeric(site_names)
J_plants[, 2] <- as.numeric(MDT)
for (spe_p in 1:no_sites)
{
  H_p <- diversity(p_sitespeciesmatrix[[spe_p]])
  J_p <- H / log(length(p_sitespeciesmatrix[[spe_p]]))
  J_plants[spe_p, 3] <- as.numeric(J_p)
}

speciesevenness <-
  cbind(as.character(site_names), MDT, J_insects[, 3], J_plants[, 3])
write.table(speciesevenness, paste(OUT_DIR, "insects-plants_sp-evenness.txt", sep = "/"), sep = "\t")

##INTERACTION EVENNESS CALCUALTIONS

#libraries
library(bipartite)
library(lme4)
library(ggplot2)
library(lattice)
library(grid)
library(gridBase)
library(plotrix)

#import and attach interaction data
allwebdata <-   read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = T, sep = ",", quote = "")
attach(allwebdata)

#site names ordered alphabetically
site_names <- sort(unique(Site))

#MDT
MDT <- unique(M_D_T)

#create interaction data matrix for each site for all bipartite networks 
#together
allinsects_datamatrix <-
  frame2webs(
    allwebdata,
    varnames = c("Lower_Taxon", "Upper_Taxon", "Site"),
    emptylist = TRUE
  )


lower_webs <-
  subset(
    allwebdata,
    allwebdata$Web == "PL_FV" |
      allwebdata$Web == "PL_LM" |
      allwebdata$Web == "PL_CP" | allwebdata$Web == "PL_SF"
  )
upper_webs <-
  subset(
    allwebdata,
    allwebdata$Web == "LM_P" |
      allwebdata$Web == "CP_P" | allwebdata$Web == "SF_P"
  )
lower_webs_datamatrix <-
  frame2webs(
    lower_webs,
    varnames = c("Lower_Taxon", "Upper_Taxon", "Site"),
    emptylist = TRUE
  )
upper_webs_datamatrix <-
  frame2webs(
    upper_webs,
    varnames = c("Lower_Taxon", "Upper_Taxon", "Site"),
    emptylist = TRUE
  )

no_sites <- length(site_names) #number of sites

IE_all_insects <-
  matrix(nrow = 30, ncol = 2) #create empty matrix for interaction evenness

for (i in 1:no_sites) {
  IE_all_insects[i, 1] <-
    as.character(site_names[i]) #put site name into matrix
  IE_all_insects[i, 2] <-
    networklevel(allinsects_datamatrix[[i]],
                 index = "interaction evenness",
                 weighted = TRUE) #calculate interaction evenness
}
write.table(IE_all_insects, file = paste(OUT_DIR, "all-webs_interaction-evenness.txt", sep = "/"), sep = "\t")

detach(allwebdata)

###MANOVA

#import and attach interaction data
site_metrics_data <- read.table(paste(DAT_DIR, "plant_abundances", "extracted_site_level_data.txt", sep = "/"), header = T, sep = "\t", quote = "")
attach(site_metrics_data)

#MANOVA
Y <-
  cbind(
    plant_rich,
    log_FU_abundance,
    total_insect_sp_rich,
    total_insect_sp_abund,
    (IE),
    (species_even_insects)
  )
fit_all <- manova(Y ~ MDT, data = site_metrics_data)
summary.manova(fit_all, test = "Pillai")
summary.aov(fit_all)

detach(site_metrics_data)

#post-hocs
md <- subset(site_metrics_data, site_metrics_data$MDT != 3)
attach(md)
Y_md <-
  cbind(
    plant_rich,
    log_FU_abundance,
    total_insect_sp_rich,
    total_insect_sp_abund,
    (IE),
    (species_even_insects)
  )
fit_md <- manova(Y_md ~ MDT, data = md)
summary(fit_md, test = "Pillai")
summary.aov(fit_md)
detach(md)

mt <- subset(site_metrics_data, site_metrics_data$MDT != 2)
attach(mt)
Y_mt <-
  cbind(
    plant_rich,
    log_FU_abundance,
    total_insect_sp_rich,
    total_insect_sp_abund,
    (IE),
    (species_even_insects)
  )
fit_mt <- manova(Y_mt ~ MDT, data = mt)
summary(fit_mt, test = "Pillai")
summary.aov(fit_mt)
detach(mt)

dt <- subset(site_metrics_data, site_metrics_data$MDT != 1)
attach(dt)
Y_dt <-
  cbind(
    plant_rich,
    log_FU_abundance,
    total_insect_sp_rich,
    total_insect_sp_abund,
    (IE),
    (species_even_insects)
  )
fit_dt <- manova(Y_dt ~ MDT, data = dt)
summary(fit_dt, test = "Pillai")
summary.aov(fit_dt)
detach(dt)


#load library
library(mcglm)

#import and attach interaction data
site_metrics_data <- read.table(paste(DAT_DIR, "plant_abundances", "extracted_site_level_data.txt", sep = "/"), header = T, sep = "\t", quote = "")
attach(site_metrics_data)

#models
#plant richness
plant_rich_model <-
  glm(plant_rich ~ MDT , data = site_metrics_data, family = "poisson")
summary(plant_rich_model)
plot(plant_rich_model)
#Floral unit abundance
FU_model <-
  glm(log_FU_abundance ~ MDT , data = site_metrics_data, family = "poisson")
summary(FU_model)
plot(FU_model)
#insect species richness
insect_rich_model <-
  glm(total_insect_sp_rich ~ MDT , data = site_metrics_data, family = "poisson")
summary(insect_rich_model)
plot(insect_rich_model)
#insect abundance
insect_abundance_model <-
  glm(total_insect_sp_abund ~ MDT ,
      data = site_metrics_data,
      family = "poisson")
summary(insect_abundance_model)
plot(insect_abundance_model)
#Interaction evenness
IE_model <- glm(IE ~ MDT , data = site_metrics_data, family = "gaussian")
summary(IE_model)
plot(IE_model)
#Species evenness
SPE_model <-
  glm((species_even_insects ^ 2) ~ MDT ,
      data = site_metrics_data,
      family = "gaussian")
summary(SPE_model)
plot(SPE_model)
