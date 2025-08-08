# Mean and standard deviation of robustness calculated for the 30 sites of the landscape project (Univ. of Bristol)
# Created by Alix SAUVE on 4th May 2022

DAT_DIR <- "../outputs"

robout_files <- list.files(DAT_DIR, pattern = "_Robustness.csv")
sites_names <- unlist(strsplit(robout_files, "_Robustness.csv"))

allsites_outputs <- data.frame(site = character(), MDT = character(), scenario = character(), extthrs = numeric(), flexthrs = numeric(),
                               meanrob = numeric(), sdrob = numeric())

rob_scenario <- data.frame(scenario = c("most to least common", "least to most common", "random"),
                          abbrev = c("RMComm", "RLComm", "RRand"))

for (site in sites_names){
  output_site <- read.csv(file = paste0(DAT_DIR, site, "_Robustness.csv"))
  
  ext_thrs <- unique(output_site$ExtThrs); flex_thrs <- unique(output_site$FlexThrs)
  
  for (ext in ext_thrs){
    for (flex in flex_thrs){
      index_extflex <- (output_site$ExtThrs == ext) & (output_site$FlexThrs == flex)
      for (s in 1:nrow(rob_scenario)){
        output_site_case <- data.frame(site = site, MDT = output_site$MDT[1], scenario = rob_scenario$scenario[s],
                                       extthrs = ext, flexthrs = flex, meanrob = NA, sdrob = NA)
        meanrob_case <- mean(output_site[index_extflex, as.character(rob_scenario$abbrev[s])])
        sdrob_case <- sd(output_site[index_extflex, as.character(rob_scenario$abbrev[s])])
        output_site_case$meanrob <- meanrob_case
        output_site_case$sdrob <- sdrob_case
        
        allsites_outputs <- rbind(allsites_outputs,
                                  output_site_case)
      }
    }
  }
}
write.table(allsites_outputs, file = paste(DAT_DIR, "allsites_robustness_meansd.csv", sep = "/"), quote = FALSE, row.names = FALSE, sep = ",")
