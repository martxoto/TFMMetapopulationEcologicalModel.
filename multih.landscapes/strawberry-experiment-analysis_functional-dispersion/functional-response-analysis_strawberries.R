#load libraries
library(lmerTest)

DAT_DIR <- "../data" # change path to your working directory if need be
OUT_DIR <- "outputs"

#Strawberry data
strawdata <- read.table(paste(DAT_DIR, "strawberry_experiment", "strawberry-weight-class_1m.txt", sep = "/"), header = T, sep = "\t")
attach(strawdata)
monads <- subset(strawdata, strawdata$MorT == "M")
triads <- subset(strawdata, strawdata$MorT == "T")

#Histograms
hist(
  monads$Fruit_weight,
  col = "red",
  xlim = c(0, 3),
  ylim = c(0, 30)
)
par(new = TRUE)
hist(
  triads$Fruit_weight,
  col = "blue",
  xlim = c(0, 3),
  ylim = c(0, 30)
)

#Mixed effect model for fruit weight with site as a random factor and number of
#habitats as fixed
fruit_weight_model <-
  lmer(Fruit_weight ~ Habitats + (1 | Site), data = strawdata)
summary(fruit_weight_model)
plot(fruit_weight_model)

#Boxplot for fruit weight
boxplot(monads$Fruit_weight,
        triads$Fruit_weight,
        ylab = "Fruit Weight (g)",
        xaxt = "n")
axis(1, at = c(1, 2), labels = c("Monads", "Triads"))

detach(strawdata)

##fruit class analysis

#read in and attach data
strawclassdata <- read.table(paste(DAT_DIR, "strawberry_experiment", "strawberry-proportion-class.txt", sep = "/"), header = T, sep = "\t")
attach(strawclassdata)
Habitats <- factor(Habitats)

#subset data
monads_st <- subset(strawclassdata, strawclassdata$Habitats == 1)
triads_st <- subset(strawclassdata, strawclassdata$Habitats == 3)

#t-test for proportion class 1
t.test(monads_st$Proportion1_1mon, triads_st$Proportion1_1mon)

#Boxplot for fruit class
boxplot(
  monads_st$Proportion1_1mon,
  triads_st$Proportion1_1mon,
  ylab = "Proportion of class 1 strawberries",
  xaxt = "n"
)
axis(1, at = c(1, 2), labels = c("Monads", "Triads"))

detach(strawclassdata)

#########################################################################
# -------------------- Functional Diversity Analysis -------------------- 

# Author: Kate Pereira Maia
# Date: 06/2021

# Part 1: Interaction matrix - defines the interaction matrix to be used.
# Part 2: Calculates distances - defines the distance matrix to be used.
# Part 3: Evaluates functional space (code from Maire et al. 2015 - CHECKED).
# Part 4: Calculates functional diversity.
# Part 5: Plots and Statistical analyses.

# ---------------------- Loading data and packages ----------------------

library(bipartite); library(FD); library(vegan) 

source("./functions/remove_abundance.R")
source("./functions/bespoke_distance.R")
source("./functions/mSD.R")

# Subset flower-Visitor interactions
FVwebdata <- subset(allwebdata, allwebdata$Web == "PL_FV") 

# Site names ordered alphabetically
site_data$Site <- gsub(" ", "_", site_data$Site); site_data$Site[21] <- "Pembrey"
site_names <- sort(site_data$Site)

# --------------------- Part 1 - Interaction Matrix ---------------------

# Combining sites to create a pooled interaction matrix (METAWEB)
FV_metaweb <- frame2webs(FVwebdata, varnames = c("Lower_Taxon", "Upper_Taxon", "Number"), emptylist = TRUE)[[1]]; dim(FV_metaweb) # dim: 150 524
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

# removing dyads to match strawberry analysis
data <- data10D[data10D$MDT %in% c("M", "T"),]
data$MDT <- as.numeric(factor(data$MDT))
summary(aov(abFDis ~ MDT,data = data))
t.test(data$MDT, data$abFDis)

dev.off()
boxplot(abFDis ~ MDT, data = data10D, xlab = "Number of Habitats", 
        ylab = "Functional Diversity", ylim = c(0.15, 0.7))
text(1:3, 0.65, c("a", "ab", "b"))


######

#combining and sorting strawberry data and functional dispersion data
average_fruit_weight <-
  aggregate(strawdata$Fruit_weight,
            by = list(strawdata$Site),
            FUN = mean)
FD_data <- subset(data10D, data10D$MDT != "D")
FD_data_2 <- FD_data[-c(10, 16), ]
straw_class <- strawclassdata[order(strawclassdata$Site), ]
straw_class_2 <- straw_class[-c(11, 1), ]
com_straw_data <-
  cbind(
    as.character(average_fruit_weight$Group.1),
    as.character(FD_data_2$MDT),
    as.numeric(average_fruit_weight$x),
    as.numeric(straw_class_2$Proportion1_1mon),
    as.numeric(FD_data_2$abFDis)
  )

#linear models to determine if Functional Dispersion predicts Fruit weight and
#Fruit class
plot(as.numeric(com_straw_data[, 3]) ~ com_straw_data[, 5])
plot(com_straw_data[, 4] ~ com_straw_data[, 5])
straw_model_fruitweight <-
  lm(as.numeric(com_straw_data[, 3]) ~ as.numeric(com_straw_data[, 5]) * 
       com_straw_data[, 2])
summary(straw_model_fruitweight)
plot(straw_model_fruitweight)

straw_model_fruitclass <-
  lm(as.numeric(com_straw_data[, 4]) ~ as.numeric(com_straw_data[, 5]) * 
       com_straw_data[, 2])
summary(straw_model_fruitclass)
plot(straw_model_fruitclass)

#calculate magnitude of pollination effect
M_T_fruit_weight <-
  aggregate(as.numeric(com_straw_data[, 3]),
            by = list(com_straw_data[, 2]),
            FUN = mean)
fruit_weight_percent <-
  ((M_T_fruit_weight[2, 2] - M_T_fruit_weight[1, 2]) / M_T_fruit_weight[2, 2])
#0.05781617 or 5.8% increase
M_T_fruit_class <-
  aggregate(as.numeric(com_straw_data[, 4]),
            by = list(com_straw_data[, 2]),
            FUN = mean)
fruit_class_percent <-
  (((M_T_fruit_class[2, 2] - M_T_fruit_class[1, 2])) / M_T_fruit_class[2, 2])
#0.303244 or 30.3%

#strawberry figures
MTcol <-
  c(
    "#91bfdb",
    "#ffffbf",
    "#91bfdb",
    "#ffffbf",
    "#91bfdb",
    "#ffffbf",
    "#ffffbf",
    "#91bfdb",
    "#91bfdb",
    "#ffffbf",
    "#91bfdb",
    "#91bfdb",
    "#91bfdb",
    "#ffffbf",
    "#ffffbf",
    "#91bfdb",
    "#91bfdb",
    "#ffffbf"
  )
MTshapes <-
  c(24, 22, 24, 22, 24, 22, 22, 24, 24, 22, 24, 24, 24, 22, 22, 24, 24, 22)

com_straw_data_M <- subset(com_straw_data, com_straw_data[,2]=="M")
com_straw_data_T <- subset(com_straw_data, com_straw_data[,2]=="T")


pdf(paste(OUT_DIR, "strawberry-combined-figures_4in1.pdf", sep = "/"))

par(mfrow = c(2, 2), mar = c(4, 4, 1, 2))

boxplot(monads$Fruit_weight, triads$Fruit_weight,
        ylab = "Fruit Weight (g)")
axis(
  side = 1,
  at = c(1, 2),
  labels = c("Monad", "Triad")
)
text(x = 0.5,
     y = 2.7,
     labels = substitute(paste(bold("A"))))

boxplot(monads_st$Proportion1_1mon,
        triads_st$Proportion1_1mon,
        ylab = "Proportion of class 1 strawberries")
axis(
  side = 1,
  at = c(1, 2),
  labels = c("Monad", "Triad")
)
text(x = 0.5,
     y = 0.98,
     labels = substitute(paste(bold("B"))))

boxplot(as.numeric(com_straw_data_M[,5]),
        as.numeric(com_straw_data_T[,5]),
        xaxt = "n",
        ylab = "Functional Dispersion")
axis(
  side = 1,
  at = c(1, 2),
  labels = c("Monad", "Triad")
)
text(x = 0.5,
     y = 0.6,
     labels = substitute(paste(bold("C"))))

plot(
  com_straw_data[, 4] ~ com_straw_data[, 5],
  xlab = "Functional Dispersion",
  ylab = "Proportion of class 1 strawberries",
  pch = MTshapes,
  bg = MTcol
)
abline(lm(as.numeric(com_straw_data[, 4]) ~ as.numeric(com_straw_data[, 5])))
text(x = 0.277,
     y = 0.98,
     labels = substitute(paste(bold("D"))))

dev.off()
