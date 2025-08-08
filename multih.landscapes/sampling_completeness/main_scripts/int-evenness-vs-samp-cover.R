# interaction evenness vs. sampling cover plot

rm(list = ls(all = TRUE))

library(iNEXT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(bipartite)
library(ggpubr)

# define working directory
DAT_DIR <- "../../data"
OUT_DIR <- "../outputs"

# load data
all_obs <- read.table(paste(DAT_DIR, "all_web_interactions.csv", sep = "/"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
sites_list <- sort(unique(all_obs$Site))
sitetypes_list <- unique(all_obs$M_D_T); sitetypes_list <- sitetypes_list[c(2, 1, 3)]
num_rep <- 20 # number of replicates per level of sampling coverage

# for each site
# rarefy the number of interations
# calculate the corresponding interaction evenness
# output table : SITE | MDT | SAMPL_COVER | REP_NUM | INT_EVE

df_rarefy <- data.frame(SITE = character(), MDT = character(), SAMPL_COVER = numeric(), NINT_SAMPL = numeric(), REP_NUM = numeric(), INT_EVE = numeric())
df_obs <- data.frame(SITE = sites_list, MDT = NA, SAMPL_COVER = NA, NINT_SAMPL = NA, INT_EVE = NA)
SC_obs <- c()
df_nint_SC60 <- data.frame(SITE = character(), MDT = character(), NINT_2SAMPL = numeric())

for (site in sites_list){

  site_type <- unique(all_obs$M_D_T[all_obs$Site == site]) ; df_obs$MDT[df_obs$SITE == site] <- site_type
  int_site <- all_obs %>% filter(Site == site) %>% group_by(Lower_Taxon, Upper_Taxon, Web) %>% count() %>% ungroup()
  mat_siteweb <- xtabs(n ~ Lower_Taxon + Upper_Taxon, data = int_site) # convert into a matrix
  inteve_siteweb <- networklevel(mat_siteweb, index = "interaction evenness", weighted = TRUE)
  df_obs$INT_EVE[df_obs$SITE == site] <- inteve_siteweb
  
  out_inc <- iNEXT(int_site$n, datatype = "abundance", knots = 50, q = 0)
  tb_rarefy_siteweb <- out_inc$iNextEst$coverage_based[out_inc$iNextEst$coverage_based$Method != "Extrapolation", ]
  tb_rarefy_siteweb <- tb_rarefy_siteweb[(tb_rarefy_siteweb$m > 1), c("SC", "m")]
  
  df_obs$SAMPL_COVER[df_obs$SITE == site] <- out_inc$DataInfo$SC
  df_obs$NINT_SAMPL[df_obs$SITE == site] <- out_inc$DataInfo$n
  # print(nrow(tb_rarefy_siteweb))
  # print(out_inc$DataInfo$SC)
  SC_obs <- c(SC_obs, out_inc$DataInfo$SC)
  # df_rarefy_site <- df_rarefy[0, ]
  # for (rep in 1:num_rep){
  #   for (n in 1:nrow(tb_rarefy_siteweb)){
  #     N_int2samp <- tb_rarefy_siteweb$m[n]
  #     int2samp <- sample(1:nrow(int_site), N_int2samp, replace = TRUE, prob = int_site$n)
  # 
  #     rarefied_siteweb <- int_site[int2samp, c("Lower_Taxon", "Upper_Taxon")]
  #     rarefied_siteweb <- rarefied_siteweb %>% group_by(Lower_Taxon, Upper_Taxon) %>% count() %>% ungroup()
  #     mat_rarefied_siteweb <- xtabs(n ~ Lower_Taxon + Upper_Taxon, data = rarefied_siteweb) # convert into a matrix
  #     inteve_rarefied_siteweb <- networklevel(mat_rarefied_siteweb, index = "interaction evenness", weighted = TRUE)
  # 
  # 
  #     df_rarefy_site <- rbind(df_rarefy_site,
  #                             data.frame(SITE = site, MDT = site_type, SAMPL_COVER = tb_rarefy_siteweb$SC[n], NINT_SAMPL = N_int2samp,
  #                                        REP_NUM = rep, INT_EVE = inteve_rarefied_siteweb))
  #   }
  # }
  # df_rarefy <- rbind(df_rarefy, df_rarefy_site)
  m_SC60 <- estimateD(int_site$n, datatype = "abundance", level = 0.6, q = 0, base = "coverage")$m
  df_nint_SC60 <- rbind(df_nint_SC60,
                          data.frame(SITE = site, MDT = site_type, NINT_2SAMP = m_SC60))
}

write.csv(df_rarefy, file = paste(OUT_DIR, "df_int-evenness_rarefied-site.csv", sep = "/"), row.names = FALSE, quote = FALSE)
write.csv(df_rarefy_SC60, file = paste(OUT_DIR, "df_nint_rarefied-site60.csv", sep = "/"), row.names = FALSE, quote = FALSE)

for (site in sites_list){
  df_rarefy_site <- df_rarefy[df_rarefy$SITE == site, ]
  test_site <- cor.test(df_rarefy_site$SAMPL_COVER, df_rarefy_site$INT_EVE)
  if (test_site$p.value <= 0.05){
    print(site)
  }
}

for (sitetype in sitetypes_list){
  df_rarefy_sitetype <- df_rarefy[df_rarefy$MDT == sitetype, ]
  df_obs_sitetype <- df_obs[df_obs$MDT == sitetype, ]
  
  p.vals <- sapply(unique(df_rarefy_sitetype$SITE), function(i) {
    coef(summary(lm(INT_EVE ~ SAMPL_COVER, data = df_rarefy_sitetype[df_rarefy_sitetype$SITE == i, ])))[2,4]
  })
  
  plot_sitetype <- ggplot(df_rarefy_sitetype, aes(x = SAMPL_COVER, y = INT_EVE, group = SITE, color = SITE)) +
    geom_point(size = 1) + xlim(0, 1) + ylim(0.3, 0.9) +
    labs(title = paste0(toupper(sitetype), "S"), x = "Sampling coverage", y = "Interaction evenness") +
    scale_color_manual(values = brewer.pal(10, "Set3")) +
    # geom_smooth(method = 'lm') +
    geom_smooth(data = df_rarefy_sitetype[df_rarefy_sitetype$SITE %in% names(p.vals)[p.vals < 0.05],], aes(SAMPL_COVER, INT_EVE, colour = SITE), method='lm') +
    # geom_hline(yintercept = df_obs_sitetype$INT_EVE, linetype = "dashed", color = brewer.pal(10, "Set3")) +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 5)) +
    guides(color = guide_legend(ncol = 3))# +
    # geom_vline(xintercept = df_obs_sitetype$SAMPL_COVER, linetype = "dashed", color = brewer.pal(10, "Set3"))
  assign(paste("plot_rarefy", sitetype, sep = "_"), plot_sitetype)
}
pdf(paste(OUT_DIR, "int-even_sampl-cover.pdf", sep = "/"), width = 12, height = 6)
ggarrange(plot_rarefy_monad, plot_rarefy_dyad, plot_rarefy_triad, ncol = 3)
dev.off()


# rarefying networks to SC = 60%
df_rarefy_SC60 <- data.frame(SITE = character(), MDT = character(), NINT_SAMPL = numeric(), REP_NUM = numeric(), INT_EVE = numeric())
for (site in sites_list){
  print(site)
  site_type <- unique(all_obs$M_D_T[all_obs$Site == site]) ; df_obs$MDT[df_obs$SITE == site] <- site_type
  int_site <- all_obs %>% filter(Site == site) %>% group_by(Lower_Taxon, Upper_Taxon, Web) %>% count() %>% ungroup()

  df_rarefy_SC60_site <- df_rarefy_SC60[0, ]
  for (rep in 1:num_rep){
    
    N_int2samp <- df_nint_SC60$NINT_2SAMP[df_nint_SC60$SITE == site]
    N_int2samp_rep <- sample(c(floor(N_int2samp), ceiling(N_int2samp)), 1)
    int2samp <- sample(1:nrow(int_site), N_int2samp_rep, replace = TRUE, prob = int_site$n)

    rarefied_siteweb <- int_site[int2samp, c("Lower_Taxon", "Upper_Taxon")]
    rarefied_siteweb <- rarefied_siteweb %>% group_by(Lower_Taxon, Upper_Taxon) %>% count() %>% ungroup()
    mat_rarefied_siteweb <- xtabs(n ~ Lower_Taxon + Upper_Taxon, data = rarefied_siteweb) # convert into a matrix
    inteve_rarefied_siteweb <- networklevel(mat_rarefied_siteweb, index = "interaction evenness", weighted = TRUE)


    df_rarefy_SC60_site <- rbind(df_rarefy_SC60_site,
                            data.frame(SITE = site, MDT = site_type, NINT_SAMPL = N_int2samp_rep,
                                       REP_NUM = rep, INT_EVE = inteve_rarefied_siteweb))
  }
  df_rarefy_SC60 <- rbind(df_rarefy_SC60, df_rarefy_SC60_site)
}

ggplot(df_rarefy_SC60, aes(x = MDT, y = INT_EVE, color = MDT)) +
  geom_boxplot() + scale_colour_viridis_d() +
  labs(y = "Interaction evenness", x = "Landscape type") + theme(legend.position =  "none")

write.csv(df_rarefy_SC60, file = paste(OUT_DIR, "df_int-evenness_rarefied-site_sc60.csv", sep = "/"))

df_rarefy_SC60 <- df_rarefy_SC60 %>% group_by(SITE, MDT) %>% summarise(MEAN_INT_EVE = mean(INT_EVE), SD_INT_EVE = sd(INT_EVE))
df_rarefy_SC60$MDT <- factor(df_rarefy_SC60$MDT, levels = c("monad", "dyad", "triad"))

pdf(paste(OUT_DIR, "int-even_sampl-cover-60.pdf", sep = "/"), width = 5, height = 6)
ggplot(df_rarefy_SC60, aes(x = MDT, y = MEAN_INT_EVE, color = MDT)) +
  geom_boxplot() + scale_colour_viridis_d() +
  labs(y = "Interaction evenness", x = "Landscape type") + theme(legend.position =  "none")
dev.off()

