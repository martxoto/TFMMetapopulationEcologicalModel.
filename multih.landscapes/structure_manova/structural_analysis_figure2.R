# +++++++++++++++++++++++++++++++++++++++++++
#### Tidied version for Nature, Mar 2024 ####

# Code for Figure 2: structural analyses

#load required libraries
library(lattice)
library(grid)
library(gridBase)
library(plotrix)

DAT_DIR <- "../data" # change path to your working directory if need be
OUT_DIR <- "outputs"

#read data
site_data <- read.table(paste(DAT_DIR, "plant_abundances", "extracted_site_level_data.txt", sep = "/"), header = T, sep = "\t", quote = "")

attach(site_data)

#habitat colours
HabList <-
  c("woodland",
    "salt_marsh",
    "sand_dune",
    "grassland",
    "scrub",
    "heathland")
ColSetAll <-
  c("#fc8d59",
    "#e0f3f8",
    "#fee090",
    "#4575b4",
    "#d73027",
    "#91bfdb")

#subset MDT
m <- subset(site_data, site_data$MDT == 1)
d <- subset(site_data, site_data$MDT == 2)
t <- subset(site_data, site_data$MDT == 3)
FU_m <- log1p(m$openFU_abundance)
FU_d <- log1p(d$openFU_abundance)
FU_t <- log1p(t$openFU_abundance)

#plot figure
pdf(paste(OUT_DIR, "descriptive_figures.pdf", sep = "/"))
par(mfrow = c(3, 2), mar = c(2, 4, 1, 2))

#floral units
boxplot(
  FU_m,
  FU_d,
  FU_t,
  ylab = "log(Number of open floral units+1)",
  xaxt = "n",
  pch = 16,
  lty = 1,
  outline = FALSE
)
for (mh in 1:length(FU_m))
{
  par(new = TRUE)
  pie_mh <- 1
  y_pie <- FU_m[mh]
  pie_col <-
    c(
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[3],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  floating.pie(
    xpos = jitter(1, factor = 6),
    y = y_pie,
    pie_mh,
    radius = .06,
    col = pie_col[mh],
    startpos = 3.14 / 2
  )
}
for (dh in 1:length(FU_d))
{
  par(new = TRUE)
  pie_dh <- c(0.5, 0.5)
  y_pie <- FU_d[dh]
  pie_col_1 <-
    c(
      ColSetAll[1],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[2]
    )
  floating.pie(
    xpos = jitter(2, factor = 3),
    y = y_pie,
    pie_dh,
    radius = .06,
    col = c(pie_col_1[dh], pie_col_2[dh]),
    startpos = 3.14 / 2
  )
}
for (th in 1:length(FU_t))
{
  par(new = TRUE)
  y_pie <- FU_t[th]
  pie_th <- c(0.3, 0.3, 0.3)
  pie_col_1 <-
    c(
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[5]
    )
  pie_col_3 <-
    c(
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[3],
      ColSetAll[5],
      ColSetAll[3]
    )
  floating.pie(
    xpos = jitter(3, factor = 2),
    y = y_pie,
    pie_th,
    radius = .06,
    col = c(pie_col_1[th], pie_col_2[th], pie_col_3[th])
  )
}
text(0.5, 11.5, "A")

#plant species richness
boxplot(
  plant_rich ~ MDT,
  ylab = "Number of plant species",
  xaxt = "n",
  xlab = "",
  pch = 16,
  lty = 1,
  outline = FALSE
)
for (mh in 1:length(m$plant_rich))
{
  par(new = TRUE)
  pie_mh <- 1
  pie_col <-
    c(
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[3],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  floating.pie(
    xpos = jitter(1, factor = 6),
    y = m$plant_rich[mh],
    pie_mh,
    radius = .06,
    col = pie_col[mh],
    startpos = 3.14 / 2
  )
}
for (dh in 1:length(d$plant_rich))
{
  par(new = TRUE)
  pie_dh <- c(0.5, 0.5)
  pie_col_1 <-
    c(
      ColSetAll[1],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[2]
    )
  floating.pie(
    xpos = jitter(2, factor = 3),
    y = d$plant_rich[dh],
    pie_dh,
    radius = .06,
    col = c(pie_col_1[dh], pie_col_2[dh]),
    startpos = 3.14 / 2
  )
}
for (th in 1:length(t$plant_rich))
{
  par(new = TRUE)
  pie_th <- c(0.3, 0.3, 0.3)
  pie_col_1 <-
    c(
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[5]
    )
  pie_col_3 <-
    c(
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[3],
      ColSetAll[5],
      ColSetAll[3]
    )
  floating.pie(
    xpos = jitter(3, factor = 2),
    y = t$plant_rich[th],
    pie_th,
    radius = .06,
    col = c(pie_col_1[th], pie_col_2[th], pie_col_3[th])
  )
}
text(0.5, 87.5, "B")

#Number of insects
boxplot(
  total_insect_sp_abund ~ MDT,
  ylab = "Number of insects",
  xaxt = "n",
  xlab = "",
  pch = 16,
  lty = 1,
  outline = FALSE
)
for (mh in 1:length(m$total_insect_sp_abund))
{
  par(new = TRUE)
  pie_mh <- 1
  pie_col <-
    c(
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[3],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  floating.pie(
    xpos = jitter(1, factor = 6),
    y = m$total_insect_sp_abund[mh],
    pie_mh,
    radius = .06,
    col = pie_col[mh],
    startpos = 3.14 / 2
  )
}
for (dh in 1:length(d$total_insect_sp_abund))
{
  par(new = TRUE)
  pie_dh <- c(0.5, 0.5)
  pie_col_1 <-
    c(
      ColSetAll[1],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[2]
    )
  floating.pie(
    xpos = jitter(2, factor = 3),
    y = d$total_insect_sp_abund[dh],
    pie_dh,
    radius = .06,
    col = c(pie_col_1[dh], pie_col_2[dh]),
    startpos = 3.14 / 2
  )
}
for (th in 1:length(t$total_insect_sp_abund))
{
  par(new = TRUE)
  pie_th <- c(0.3, 0.3, 0.3)
  pie_col_1 <-
    c(
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[5]
    )
  pie_col_3 <-
    c(
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[3],
      ColSetAll[5],
      ColSetAll[3]
    )
  floating.pie(
    xpos = jitter(3, factor = 2),
    y = t$total_insect_sp_abund[th],
    pie_th,
    radius = .06,
    col = c(pie_col_1[th], pie_col_2[th], pie_col_3[th])
  )
}
text(0.5, 1075, "C")

#insect species richness
boxplot(
  total_insect_sp_rich ~ MDT,
  ylab = "Number of insect species",
  xaxt = "n",
  xlab = "",
  pch = 16,
  lty = 1,
  outline = FALSE
)
for (mh in 1:length(m$total_insect_sp_rich))
{
  par(new = TRUE)
  pie_mh <- 1
  pie_col <-
    c(
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[3],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  floating.pie(
    xpos = jitter(1, factor = 6),
    y = m$total_insect_sp_rich[mh],
    pie_mh,
    radius = .06,
    col = pie_col[mh],
    startpos = 3.14 / 2
  )
}
for (dh in 1:length(d$total_insect_sp_rich))
{
  par(new = TRUE)
  pie_dh <- c(0.5, 0.5)
  pie_col_1 <-
    c(
      ColSetAll[1],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[2]
    )
  floating.pie(
    xpos = jitter(2, factor = 3),
    y = d$total_insect_sp_rich[dh],
    pie_dh,
    radius = .06,
    col = c(pie_col_1[dh], pie_col_2[dh]),
    startpos = 3.14 / 2
  )
}
for (th in 1:length(t$total_insect_sp_rich))
{
  par(new = TRUE)
  pie_th <- c(0.3, 0.3, 0.3)
  pie_col_1 <-
    c(
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[5]
    )
  pie_col_3 <-
    c(
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[3],
      ColSetAll[5],
      ColSetAll[3]
    )
  floating.pie(
    xpos = jitter(3, factor = 2),
    y = t$total_insect_sp_rich[th],
    pie_th,
    radius = .06,
    col = c(pie_col_1[th], pie_col_2[th], pie_col_3[th])
  )
}
text(0.5, 185, "D")

#Interaction evenness
boxplot(
  IE ~ MDT,
  ylab = "Interaction evenness",
  xaxt = "n",
  pch = 16,
  lty = 1,
  outline = FALSE
)
axis(1,
     at = c(1, 2, 3),
     labels = c("Monads", "Dyads", "Triads"))
for (mh in 1:length(m$IE))
{
  par(new = TRUE)
  pie_mh <- 1
  pie_col <-
    c(
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[3],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  floating.pie(
    xpos = jitter(1, factor = 6),
    y = m$IE[mh],
    pie_mh,
    radius = .06,
    col = pie_col[mh],
    startpos = 3.14 / 2
  )
}
for (dh in 1:length(d$IE))
{
  par(new = TRUE)
  pie_dh <- c(0.5, 0.5)
  pie_col_1 <-
    c(
      ColSetAll[1],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[2]
    )
  floating.pie(
    xpos = jitter(2, factor = 3),
    y = d$IE[dh],
    pie_dh,
    radius = .06,
    col = c(pie_col_1[dh], pie_col_2[dh]),
    startpos = 3.14 / 2
  )
}
for (th in 1:length(t$IE))
{
  par(new = TRUE)
  pie_th <- c(0.3, 0.3, 0.3)
  pie_col_1 <-
    c(
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[5]
    )
  pie_col_3 <-
    c(
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[3],
      ColSetAll[5],
      ColSetAll[3]
    )
  floating.pie(
    xpos = jitter(3, factor = 2),
    y = t$IE[th],
    pie_th,
    radius = .06,
    col = c(pie_col_1[th], pie_col_2[th], pie_col_3[th])
  )
}
text(0.5, 0.6075, "E")

#species evenness
boxplot(
  species_even_insects ~ MDT,
  ylab = "Species evenness",
  xaxt = "n",
  pch = 16,
  lty = 1,
  outline = FALSE
)
axis(1,
     at = c(1, 2, 3),
     labels = c("Monads", "Dyads", "Triads"))
for (mh in 1:length(m$species_even_insects))
{
  par(new = TRUE)
  pie_mh <- 1
  pie_col <-
    c(
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[3],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  floating.pie(
    xpos = jitter(1, factor = 6),
    y = m$species_even_insects[mh],
    pie_mh,
    radius = .06,
    col = pie_col[mh],
    startpos = 3.14 / 2
  )
}
for (dh in 1:length(d$species_even_insects))
{
  par(new = TRUE)
  pie_dh <- c(0.5, 0.5)
  pie_col_1 <-
    c(
      ColSetAll[1],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[3],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[2]
    )
  floating.pie(
    xpos = jitter(2, factor = 3),
    y = d$species_even_insects[dh],
    pie_dh,
    radius = .06,
    col = c(pie_col_1[dh], pie_col_2[dh]),
    startpos = 3.14 / 2
  )
}
for (th in 1:length(t$species_even_insects))
{
  par(new = TRUE)
  pie_th <- c(0.3, 0.3, 0.3)
  pie_col_1 <-
    c(
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[4],
      ColSetAll[6],
      ColSetAll[4]
    )
  pie_col_2 <-
    c(
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[1],
      ColSetAll[1],
      ColSetAll[6],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[5],
      ColSetAll[1],
      ColSetAll[5]
    )
  pie_col_3 <-
    c(
      ColSetAll[2],
      ColSetAll[1],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[5],
      ColSetAll[2],
      ColSetAll[3],
      ColSetAll[3],
      ColSetAll[5],
      ColSetAll[3]
    )
  floating.pie(
    xpos = jitter(3, factor = 2),
    y = t$species_even_insects[th],
    pie_th,
    radius = .06,
    col = c(pie_col_1[th], pie_col_2[th], pie_col_3[th])
  )
}
legend (
  "bottomright",
  legend = c(
    "Woodland",
    "Salt marsh",
    "Sand dune",
    "Grassland",
    "Scrub",
    "Heathland"
  ),
  fill = ColSetAll,
  ncol = 2
)
text(0.5, 0.89, "F")
dev.off()

detach(site_data)