# this script draws the cumulative secondary extinction curves following plants removal in each sites sampled by Hackett et al.
# edited 23/10/2021 by Alix SAUVE

rm(list=ls(all=TRUE))

# loading libraries
library(tidyr)
library(AICcmodavg)
library(car)
library(gtools)
library(lme4) # using the same package as for the test on the pollination experiment
library(lmerTest)

SIM_DIR <- "../outputs"

# get the site list
SiteList <- list.files(SIM_DIR, pattern = "csv") # list of csv files (outputs of the robustness analysis)
SiteList <- unlist(strsplit(SiteList, split = "_ExtSeq"))
SiteList <- unique(SiteList[-grep(".csv", SiteList)])

# scenario list
ScenarioList <- c("Random", "LeastComm", "MostComm")
ScenarioCol <- RColorBrewer::brewer.pal(3, "Set1"); names(ScenarioCol) <- ScenarioList

MDT <- data.frame(Site = SiteList, MDT = NA)
# for each site, draw the robustness curves
for (Site in SiteList){
  # Site = SiteList[1]
  setwd(SIM_DIR)
  # random extinction scenario
  # get the output file
  FileName <- paste(Site, "ExtSeq_Random.csv", sep = "_")
  RandomSeq <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  # least to most common plants extinction scenario
  # get the output file
  FileName <- paste(Site, "ExtSeq_LeastComm.csv", sep = "_")
  LeastCommSeq <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  # most to least common plants extinction scenario
  # get the output file
  FileName <- paste(Site, "ExtSeq_MostComm.csv", sep = "_")
  MostCommSeq <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  # robustness
  # get the output file
  FileName <- paste(Site, "Robustness.csv", sep = "_")
  Robustness <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  MDT$MDT[MDT$Site == Site] <- unique(Robustness$MDT)
  
  NPlants <- ncol(MostCommSeq) - 3 - 1 # number of plants
  PropExtPl <- seq(0, NPlants) / NPlants # proportion of extinct plants
  
  # create an empty dataframe to fill with the mean extinction sequence
  MeanExtSeq <- data.frame(Scenario = character(), Ext = numeric(), Flex = numeric(), NPlants = numeric())
  for (N in 0:NPlants){
    MeanExtSeq[, paste("Prop", N, sep = "")] <- numeric()
  }
  
  NRep <- max(unique(RandomSeq$SimNo)) # number of replicates
  FlexLevels <- unique(RandomSeq$FlexThrs); ExtLevels <- unique(RandomSeq$ExtThrs) # extinction and flexibility levels
  for (Ext in ExtLevels){
    # Ext <- ExtLevels[1]
    setwd(OUT_DIR)
    pdf(paste(Site, "_Rob_E", Ext, ".pdf", sep = ""), width = 8, height = 8)
    par(mfrow = c(length(FlexLevels), 3), mar = c(4, 4, 1, 1))
    for (Flex in FlexLevels){
      # Flex <- FlexLevels[1]
      # plot each simulated extinction sequence
      for (Scenario in ScenarioList){
        # Scenario <- ScenarioList[1]
        ExtSeq <- get(paste(Scenario, "Seq", sep = ""))
        ExtSeq_ExtFlex <- subset(ExtSeq, (FlexThrs == Flex) & (ExtThrs == Ext))

        plot(0, 0, type = "n", xlab = "% plant removed", ylab = "% remaining species", xlim = c(0, 1), ylim = c(0, 1))
        for (Sim in 1:NRep){
          PropRemainingSp <- ExtSeq_ExtFlex[Sim, -c(1:3)]
          lines(PropExtPl, PropRemainingSp, col = adjustcolor("grey", alpha.f = 0.1))
        }
        MeanExtSeq_FlexExt <- apply(ExtSeq_ExtFlex[, -c(1:3)], 2, mean)
        lines(PropExtPl, MeanExtSeq_FlexExt, lwd = 3)
        if (Flex == FlexLevels[1]){
          mtext(Scenario, 3, line = 0, font = 2)
        }
        NewLine_MeanExtSeq <- data.frame(Scenario = Scenario, Ext = Ext, Flex = Flex, NPlants = NPlants)
        for (N in 0:NPlants){
          NewLine_MeanExtSeq[, paste("Prop", N, sep = "")] <- MeanExtSeq_FlexExt[N + 1]
        }
        MeanExtSeq <- rbind(MeanExtSeq, NewLine_MeanExtSeq)
      }
    }
    dev.off()
  }
  assign(paste(Site, "MeanExtSeq", sep = "_"), MeanExtSeq)
}

# plotting robustness grouping monads, dyads and triads together (mean value)
SiteCol <- RColorBrewer::brewer.pal(10, "Set3")
MDTTypes <- c("monad", "dyad", "triad")
for (Ext in ExtLevels){
  # Ext <- ExtLevels[1]
  setwd(OUT_DIR)
  pdf(paste("AllRob_E", Ext, ".pdf", sep = ""), width = 8, height = 8)
  par(mfrow = c(length(FlexLevels), 3), mar = c(4, 4, 1, 1))
  for (Flex in FlexLevels){
    # Flex <- FlexLevels[1]
    for (Type in MDTTypes){
      # Type <- MDTTypes[1]
      Landscapes <- MDT[MDT$MDT == Type, ]
      plot(0, 0, type = "n", xlab = "% plant removed", ylab = "% remaining species", xlim = c(0, 1), ylim = c(0, 1))
      NSite <- 1
      for (Site in Landscapes$Site){
        # Site <- Landscapes$Site[1]
        # least to most common plants extinction scenario
        # get the output file
        setwd(SIM_DIR)
        FileName <- paste(Site, "ExtSeq_LeastComm.csv", sep = "_")
        LeastCommSeq <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
        LeastCommSeq <- subset(LeastCommSeq, (FlexThrs == Flex) & (ExtThrs == Ext))
        NPlants <- ncol(LeastCommSeq) - 3 - 1 # number of plants
        PropExtPl <- seq(0, NPlants) / NPlants # proportion of extinct plants
        print(c(FileName, NPlants))
        for (Sim in 1:NRep){
          PropRemainingSp <- LeastCommSeq[Sim, -c(1:3)]
          lines(PropExtPl, PropRemainingSp, col = adjustcolor(SiteCol[NSite], alpha.f = 0.01))#, lwd = 0.1)#0.1), lwd = 0.1)
        }
        NSite <- 1 + NSite
      }
      # plot the average
      NSite <- 1
      for (Site in Landscapes$Site){
        ExtSeq <- get(paste(Site, "MeanExtSeq", sep = "_"))
        ExtSeq <- ExtSeq[(ExtSeq$Scenario == "LeastComm") & (ExtSeq$Ext == Ext) & (ExtSeq$Flex == Flex), ]
        NPlants <- unique(ExtSeq$NPlants)
        PropExtPl <- seq(0, NPlants) / NPlants # proportion of extinct plants
        
        lines(PropExtPl, ExtSeq[, -c(1:4)], col = SiteCol[NSite], lwd = 3)
        NSite <- 1 + NSite
      }
    }
  }
  dev.off()
}

# plotting robustness grouping monads, dyads and triads together (mean value)
SiteCol <- RColorBrewer::brewer.pal(10, "Set3")
MDTTypes <- c("monad", "dyad", "triad")
for (Ext in ExtLevels){
  setwd(OUT_DIR)
  pdf(paste("MeanRob_E", Ext, ".pdf", sep = ""), width = 8, height = 8)
  par(mfrow = c(length(FlexLevels), 3), mar = c(4, 4, 1, 1))
  for (Flex in FlexLevels){
    for (Type in MDTTypes){
      Landscapes <- MDT[MDT$MDT == Type, ]
      plot(0, 0, type = "n", xlab = "% plant removed", ylab = "% remaining species", xlim = c(0, 1), ylim = c(0, 1))
      NSite <- 1
      for (Site in Landscapes$Site){
        ExtSeq <- get(paste(Site, "MeanExtSeq", sep = "_"))
        ExtSeq <- ExtSeq[(ExtSeq$Scenario == "LeastComm") & (ExtSeq$Ext == Ext) & (ExtSeq$Flex == Flex), ]
        NPlants <- unique(ExtSeq$NPlants)
        PropExtPl <- seq(0, NPlants) / NPlants # proportion of extinct plants
        
        lines(PropExtPl, ExtSeq[, -c(1:4)], col = SiteCol[NSite], lwd = 3)
        NSite <- 1 + NSite
      }
    }
  }
  dev.off()
}

pdf(paste(OUT_DIR, "SelectRob.pdf", sep = "/"), width = 8, height = 4)
par(mfrow = c(2, 3), mar = c(4, 4, 1, 1)); SubPlot <- 1
Ext <- 0.5
for (Flex in c(0.25, 1)){
  for (Type in MDTTypes){
    Landscapes <- MDT[MDT$MDT == Type, ]
    plot(0, 0, type = "n", xlab = "% plant removed", ylab = "% remaining species", xlim = c(0, 1), ylim = c(0, 1))
    NSite <- 1
    for (Site in Landscapes$Site){
      # least to most common plants extinction scenario
      # get the output file
      setwd(SIM_DIR)
      FileName <- paste(Site, "ExtSeq_LeastComm.csv", sep = "_")
      LeastCommSeq <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
      LeastCommSeq <- subset(LeastCommSeq, (FlexThrs == Flex) & (ExtThrs == Ext))
      NPlants <- ncol(LeastCommSeq) - 3 - 1 # number of plants
      PropExtPl <- seq(0, NPlants) / NPlants # proportion of extinct plants
      print(c(FileName, NPlants))
      for (Sim in 1:NRep){
        PropRemainingSp <- LeastCommSeq[Sim, -c(1:3)]
        lines(PropExtPl, PropRemainingSp, col = adjustcolor(SiteCol[NSite], alpha.f = 0.01))#, lwd = 0.1)#0.1), lwd = 0.1)
      }
      NSite <- 1 + NSite
    }
    # plot the average
    NSite <- 1
    for (Site in Landscapes$Site){
      ExtSeq <- get(paste(Site, "MeanExtSeq", sep = "_"))
      ExtSeq <- ExtSeq[(ExtSeq$Scenario == "LeastComm") & (ExtSeq$Ext == Ext) & (ExtSeq$Flex == Flex), ]
      NPlants <- unique(ExtSeq$NPlants)
      PropExtPl <- seq(0, NPlants) / NPlants # proportion of extinct plants
      
      lines(PropExtPl, ExtSeq[, -c(1:4)], col = SiteCol[NSite], lwd = 3)
      NSite <- 1 + NSite
    }
    legend("topleft", legend = toupper(letters[SubPlot]), cex = 1.2, bty = "n"); SubPlot <- SubPlot + 1
  }
}
dev.off()

 # for each site, draw the robustness curves, this time combining scenarios
for (Site in SiteList){
  setwd(SIM_DIR)
  # random extinction scenario
  # get the output file
  FileName <- paste(Site, "ExtSeq_Random.csv", sep = "_")
  RandomSeq <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  # least to most common plants extinction scenario
  # get the output file
  FileName <- paste(Site, "ExtSeq_LeastComm.csv", sep = "_")
  LeastCommSeq <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  # most to least common plants extinction scenario
  # get the output file
  FileName <- paste(Site, "ExtSeq_MostComm.csv", sep = "_")
  MostCommSeq <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  # robustness
  # get the output file
  FileName <- paste(Site, "Robustness.csv", sep = "_")
  Robustness <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  MDT$MDT[MDT$Site == Site] <- unique(Robustness$MDT)
  
  NPlants <- ncol(MostCommSeq) - 3 - 1 # number of plants
  PropExtPl <- seq(0, NPlants) / NPlants # proportion of extinct plants
  
  # create an empty dataframe to fill with the mean extinction sequence
  MeanExtSeq <- data.frame(Scenario = character(), Ext = numeric(), Flex = numeric(), NPlants = numeric())
  for (N in 0:NPlants){
    MeanExtSeq[, paste("Prop", N, sep = "")] <- numeric()
  }
  
  NRep <- max(unique(RandomSeq$SimNo)) # number of replicates
  FlexLevels <- unique(RandomSeq$FlexThrs); ExtLevels <- unique(RandomSeq$ExtThrs) # extinction and flexibility levels
  setwd(OUT_DIR)
  pdf(paste(Site, "_Rob.pdf", sep = ""), width = 8, height = 10)
  par(mfrow = c(length(FlexLevels), length(ExtLevels)), mar = c(4, 4, 1, 1)); SubPlot <- 1
  for (Flex in FlexLevels){
    for (Ext in ExtLevels){
      # plot each simulated extinction sequence
      plot(0, 0, type = "n", xlab = "% plant removed", ylab = "% remaining species", xlim = c(0, 1), ylim = c(0, 1))
      
      for (Scenario in ScenarioList){
        # Scenario <- ScenarioList[1]
        ExtSeq <- get(paste(Scenario, "Seq", sep = ""))
        ExtSeq_ExtFlex <- subset(ExtSeq, (FlexThrs == Flex) & (ExtThrs == Ext))
        
        for (Sim in 1:NRep){
          PropRemainingSp <- ExtSeq_ExtFlex[Sim, -c(1:3)]
          lines(PropExtPl, PropRemainingSp, col = adjustcolor(ScenarioCol[Scenario], alpha.f = 0.01), lwd = 0.1)
        }
        MeanExtSeq_FlexExt <- apply(ExtSeq_ExtFlex[, -c(1:3)], 2, mean)
        assign(paste(Site, "MeanExtSeq_FlexExt", Scenario, sep = "_"), MeanExtSeq_FlexExt)
        
        # if (Flex == FlexLevels[1]){
        #   mtext(Scenario, 3, line = 0, font = 2)
        # }
        # NewLine_MeanExtSeq <- data.frame(Scenario = Scenario, Ext = Ext, Flex = Flex, NPlants = NPlants)
        # for (N in 0:NPlants){
        #   NewLine_MeanExtSeq[, paste("Prop", N, sep = "")] <- MeanExtSeq_FlexExt[N + 1]
        # }
        # MeanExtSeq <- rbind(MeanExtSeq, NewLine_MeanExtSeq)
      }
      for (Scenario in ScenarioList){
        MeanExtSeq_FlexExt <- get(paste(Site, "MeanExtSeq_FlexExt", Scenario, sep = "_"))
        lines(PropExtPl, MeanExtSeq_FlexExt, lwd = 3, col = ScenarioCol[Scenario])
      }
      legend("topleft", legend = toupper(letters[SubPlot]), bty = "n", cex = 1.5); SubPlot <- SubPlot + 1
    }
  }
  dev.off()
}

MDT <- data.frame(Site = SiteList, MDT = NA)
# for each site, identify whether M/D/T
for (Site in SiteList){
  setwd(SIM_DIR)
  # robustness
  # get the output file
  FileName <- paste(Site, "Robustness.csv", sep = "_")
  Robustness <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
  MDT$MDT[MDT$Site == Site] <- unique(Robustness$MDT)
}

MDT$SiteName <- unlist(lapply(strsplit(as.character(MDT$Site), "_"), paste, collapse = " "))
MDT$HabNo <- 1; MDT$HabNo[MDT$MDT == "dyad"] <- 2; MDT$HabNo[MDT$MDT == "triad"] <- 3
NRep <- max(unique(Robustness$SimNo)) # number of replicates
SiteCol <- RColorBrewer::brewer.pal(10, "Set3")
MDTTypes <- c("monad", "dyad", "triad")


pdf(paste(SIM_DIR, "RobustnessComp.pdf", sep = "/"), width = 8, height = 8)
par(mfrow = c(2, 1), mar = c(6, 4, 1, 1)); SubPlot <- 1
Ext <- 0.5
for (Flex in c(0.25, 1)){
  plot(0, 0, type = "n", xlim = c(0, 35), ylim = c(0.2, 1), xlab = NA, ylab = "Robustness", xaxt = "n"); NSite <- 1
  axis(1, at = c(1:10, 13:22, 25:34), labels = NA)
  text(c(1:10, 13:22, 25:34), rep(0.1, 30), srt = 90, adj = 1, labels = MDT$SiteName[order(MDT$HabNo)], xpd = TRUE, cex = 0.7, srt = 60)
  text(c(5.5, 17.5, 29.5), rep(0.25, 3), labels = toupper(paste(MDTTypes, "s", sep = "")))
  RobustnessMDT <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
  for (Type in MDTTypes){
    Landscapes <- MDT[MDT$MDT == Type, ]
    # plot robustness
    RobustnessMDTType <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
    NCol <- 1
    for (Site in Landscapes$Site){
      # get the output file
      setwd(SIM_DIR)
      FileName <- paste(Site, "Robustness.csv", sep = "_")
      RobustnessSite <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
      RobustnessSite <- subset(RobustnessSite, (ExtThrs == Ext) & (FlexThrs == Flex))
      points(rep(NSite, NRep), RobustnessSite$RLComm, pch = 20, col = adjustcolor(SiteCol[NCol], alpha.f = 0.5))
      NSite <- NSite + 1; NCol <- NCol + 1
      RobustnessMDTType <- rbind(RobustnessMDTType, RobustnessSite)
    }
    boxplot(RobustnessMDTType$RLComm, col = NA, add = TRUE, at = NSite - 11 / 2, outline = FALSE, boxwex = 20, yaxt = "n")
    NSite <- NSite + 2
    RobustnessMDT <- rbind(RobustnessMDT, RobustnessMDTType)
  }
  legend("topleft", legend = toupper(letters[SubPlot]), bty = "n", cex = 1.5); SubPlot <- SubPlot + 1
  RobustnessMDT$MDT <- as.factor(RobustnessMDT$MDT)
  print(Flex)
}
dev.off()

NRep <- max(unique(Robustness$SimNo)) # number of replicates
SiteCol <- RColorBrewer::brewer.pal(10, "Set3")
MDTTypes <- c("monad", "dyad", "triad")
MDT$SiteName <- unlist(lapply(strsplit(as.character(MDT$Site), "_"), paste, collapse = " "))
MDT$HabNo <- 1; MDT$HabNo[MDT$MDT == "dyad"] <- 2; MDT$HabNo[MDT$MDT == "triad"] <- 3

pdf(paste(SIM_DIR, "fig3_robustness.pdf", sep = "/"), width = 8, height = 6)
layout_matrix <- matrix(c(1, 1, 1, 2, 3, 4), byrow = TRUE, nrow = 2, ncol = 3)
layout(mat = layout_matrix)
par(mar = c(6, 4, 1, 1)); SubPlot <- 1
Ext <- 0.5 ; Flex <- 1
plot(0, 0, type = "n", xlim = c(0, 35), ylim = c(0.2, 1), xlab = NA, ylab = "Robustness", xaxt = "n"); NSite <- 1
axis(1, at = c(1:10, 13:22, 25:34), labels = NA)
text(c(1:10, 13:22, 25:34), rep(0.1, 30), srt = 90, adj = 1, labels = MDT$SiteName[order(MDT$HabNo)], xpd = TRUE, cex = 0.7, srt = 60)
text(c(5.5, 17.5, 29.5), rep(0.25, 3), labels = toupper(paste(MDTTypes, "s", sep = "")))
RobustnessMDT <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
for (Type in MDTTypes){
  Landscapes <- MDT[MDT$MDT == Type, ]
  # plot robustness
  RobustnessMDTType <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
  NCol <- 1
  for (Site in Landscapes$Site){
    # get the output file
    setwd(SIM_DIR)
    FileName <- paste(Site, "Robustness.csv", sep = "_")
    RobustnessSite <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
    RobustnessSite <- subset(RobustnessSite, (ExtThrs == Ext) & (FlexThrs == Flex))
    points(rep(NSite, NRep), RobustnessSite$RLComm, pch = 20, col = adjustcolor(SiteCol[NCol], alpha.f = 0.5))
    NSite <- NSite + 1; NCol <- NCol + 1
    RobustnessMDTType <- rbind(RobustnessMDTType, RobustnessSite)
  }
  boxplot(RobustnessMDTType$RLComm, col = NA, add = TRUE, at = NSite - 11 / 2, outline = FALSE, boxwex = 20, yaxt = "n")
  NSite <- NSite + 2
  RobustnessMDT <- rbind(RobustnessMDT, RobustnessMDTType)
}
legend("topleft", legend = toupper(letters[SubPlot]), bty = "n", cex = 1.5); SubPlot <- SubPlot + 1
RobustnessMDT$MDT <- as.factor(RobustnessMDT$MDT)

leveneTest(RLComm ~ MDT, data = RobustnessMDT)
for (Type in MDTTypes){
  print(c(Flex, Type, IQR(RobustnessMDT$RLComm[RobustnessMDT$MDT == Type])))
}

for (Ext in c(0.25, 0.5, 0.75)){
  BFTest <- data.frame(FlexThrs = c(0, 0.25, 0.5, 1), Fstat_LComm = numeric(4), pValue_LComm = NA,
                       Fstat_Rand = numeric(4), pValue_Rand = NA)
  for (Flex in c(0, 0.25, 0.5, 1)){
    RobustnessMDT <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
    for (Type in MDTTypes){
      Landscapes <- MDT[MDT$MDT == Type, ]
      # plot robustness
      RobustnessMDTType <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
      NCol <- 1
      for (Site in Landscapes$Site){
        # get the output file
        setwd(SIM_DIR)
        FileName <- paste(Site, "Robustness.csv", sep = "_")
        RobustnessSite <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
        RobustnessSite <- subset(RobustnessSite, (ExtThrs == Ext) & (FlexThrs == Flex))
        RobustnessMDTType <- rbind(RobustnessMDTType, RobustnessSite)
      }
      RobustnessMDT <- rbind(RobustnessMDT, RobustnessMDTType)
    }
    BFTestExtFlex_LComm <- leveneTest(RLComm ~ MDT, data = RobustnessMDT)
    BFTest$Fstat_LComm[BFTest$FlexThrs == Flex] <- BFTestExtFlex_LComm$`F value`[1]
    BFTest$pValue_LComm[BFTest$FlexThrs == Flex] <- stars.pval(BFTestExtFlex_LComm$`Pr(>F)`[1])[1]
    BFTestExtFlex_Rand <- leveneTest(RRand ~ MDT, data = RobustnessMDT)
    BFTest$Fstat_Rand[BFTest$FlexThrs == Flex] <-  BFTestExtFlex_Rand$`F value`[1]
    BFTest$pValue_Rand[BFTest$FlexThrs == Flex] <- stars.pval(BFTestExtFlex_Rand$`Pr(>F)`[1])[1]
  }
  plot(BFTest$Fstat_LComm ~ BFTest$FlexThrs, ylim = c(1, 1500), xlab = "Diet flexibility", ylab = "Brown–Forsythe test statistic", pch = 15)
  lines(BFTest$Fstat_LComm ~ BFTest$FlexThrs) ; text(BFTest$FlexThrs, BFTest$Fstat_LComm + 50 , labels = BFTest$pValue_LComm)
  points(BFTest$Fstat_Rand ~ BFTest$FlexThrs, pch = 20) ; lines(BFTest$Fstat_Rand ~ BFTest$FlexThrs, lty = 2)
  text(BFTest$FlexThrs, BFTest$Fstat_Rand + 50 , labels = BFTest$pValue_Rand)
  legend("topleft", legend = toupper(letters[SubPlot]), bty = "n") ; SubPlot <- SubPlot + 1
}
dev.off()

# Does the number of habitats (predictive variable = MDT = fixed effect) affect community robustness (response variable = RLComm)?
# Given that 500 extinctions sequences are run for each sites, we need a mixed effect model with site as a random effect to avoid pseudoreplication.

RobustnessMDT$Site <- as.factor(RobustnessMDT$Site)
RobustnessMDT$MDT <- factor(RobustnessMDT$MDT, levels = c("monad", "dyad", "triad"))

lmeRobustness <- lmer(RLComm ~ MDT + (1 | Site), data = RobustnessMDT)
# in Fig. 3A, FlexThrs = 100%, ExtThres = 0.5
RobustnessMDT_fig3A <- subset(RobustnessMDT, (ExtThrs == 1) & (FlexThrs == 0.5))
lmeRobustness_fig3A <- lmer(RLComm ~ MDT + (1 | Site), data = RobustnessMDT)
# "There was no difference in mean robustness among monads, dyads and triads (...)."

lmeRobustness_full <- lmerTest::lmer(RLComm ~ MDT + (1 | Site), data = RobustnessMDT)
lmeRobustness_null <- lmerTest::lmer(RLComm ~ (1 | Site), data = RobustnessMDT)
anova(lmeRobustness_null, lmeRobustness_full)
anova(lmeRobustness_full)

pdf(paste(OUT_DIR, "fig11_robustness_most2leastcomm-scenario.pdf", sep = "/"), width = 11.7, height = 6.3)
par(mfrow = c(3, 4), mar = c(6, 4, 1, 1)); SubPlot <- 1
Subplot <- 1
for (Ext in c(0.25, 0.5, 0.75)){
  for (Flex in c(0, 0.25, 0.5, 1)){
    plot(0, 0, type = "n", xlim = c(0, 35), ylim = c(0.2, 1), xlab = NA, ylab = "Robustness", xaxt = "n"); NSite <- 1
    axis(1, at = c(1:10, 13:22, 25:34), labels = NA)
    text(c(1:10, 13:22, 25:34), rep(0.1, 30), srt = 90, adj = 1, labels = MDT$SiteName[order(MDT$HabNo)], xpd = TRUE, cex = 0.7, srt = 60)
    text(c(5.5, 17.5, 29.5), rep(0.25, 3), labels = toupper(paste(MDTTypes, "s", sep = "")))
    RobustnessMDT <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
    for (Type in MDTTypes){
      Landscapes <- MDT[MDT$MDT == Type, ]
      # plot robustness
      RobustnessMDTType <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
      NCol <- 1
      for (Site in Landscapes$Site){
        # get the output file
        setwd(SIM_DIR)
        FileName <- paste(Site, "Robustness.csv", sep = "_")
        RobustnessSite <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
        RobustnessSite <- subset(RobustnessSite, (ExtThrs == Ext) & (FlexThrs == Flex))
        points(rep(NSite, NRep), RobustnessSite$RMComm, pch = 20, col = adjustcolor(SiteCol[NCol], alpha.f = 0.5))
        NSite <- NSite + 1; NCol <- NCol + 1
        RobustnessMDTType <- rbind(RobustnessMDTType, RobustnessSite)
      }
      boxplot(RobustnessMDTType$RMComm, col = NA, add = TRUE, at = NSite - 11 / 2, outline = FALSE, boxwex = 20, yaxt = "n")
      NSite <- NSite + 2
      RobustnessMDT <- rbind(RobustnessMDT, RobustnessMDTType)
    }
    legend("topleft", legend = toupper(letters[Subplot]), bty = "n")
    Subplot <- Subplot + 1
  }
}
dev.off()

pdf(paste(OUT_DIR, "fig12_robustness_random-scenario.pdf", sep = "/"), width = 11.7, height = 6.3)
par(mfrow = c(3, 4), mar = c(6, 4, 1, 1)); SubPlot <- 1
Subplot <- 1
for (Ext in c(0.25, 0.5, 0.75)){
  for (Flex in c(0, 0.25, 0.5, 1)){
    plot(0, 0, type = "n", xlim = c(0, 35), ylim = c(0.2, 1), xlab = NA, ylab = "Robustness", xaxt = "n"); NSite <- 1
    axis(1, at = c(1:10, 13:22, 25:34), labels = NA)
    text(c(1:10, 13:22, 25:34), rep(0.1, 30), srt = 90, adj = 1, labels = MDT$SiteName[order(MDT$HabNo)], xpd = TRUE, cex = 0.7, srt = 60)
    text(c(5.5, 17.5, 29.5), rep(0.25, 3), labels = toupper(paste(MDTTypes, "s", sep = "")))
    RobustnessMDT <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
    for (Type in MDTTypes){
      Landscapes <- MDT[MDT$MDT == Type, ]
      # plot robustness
      RobustnessMDTType <- data.frame(Site = character(), SimNo = numeric(), RMComm = numeric(), RLComm = numeric(), RRand = numeric()) # robustness data frame
      NCol <- 1
      for (Site in Landscapes$Site){
        # get the output file
        setwd(SIM_DIR)
        FileName <- paste(Site, "Robustness.csv", sep = "_")
        RobustnessSite <- read.csv(FileName, header = TRUE, stringsAsFactors = FALSE)
        RobustnessSite <- subset(RobustnessSite, (ExtThrs == Ext) & (FlexThrs == Flex))
        points(rep(NSite, NRep), RobustnessSite$RRand, pch = 20, col = adjustcolor(SiteCol[NCol], alpha.f = 0.5))
        NSite <- NSite + 1; NCol <- NCol + 1
        RobustnessMDTType <- rbind(RobustnessMDTType, RobustnessSite)
      }
      boxplot(RobustnessMDTType$RRand, col = NA, add = TRUE, at = NSite - 11 / 2, outline = FALSE, boxwex = 20, yaxt = "n")
      NSite <- NSite + 2
      RobustnessMDT <- rbind(RobustnessMDT, RobustnessMDTType)
    }
    legend("topleft", legend = toupper(letters[Subplot]), bty = "n")
    Subplot <- Subplot + 1
  }
}
dev.off()


