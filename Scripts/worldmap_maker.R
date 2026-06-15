library(ggplot2) #ggplot2_3.5.2
library(scatterpie) #scatterpie_0.2.6
library(tidyr) #tidyr_1.3.1
library(readr) #readr_2.1.5

aqui_coord <- read_delim("C:/Users/joaoa/Ambiente de Trabalho/Work/Current/Paper_Aquimarinas/Paper_Aq135/WorldMap/map_to_csv_m_vs_g.csv",delim="\t")
aqui_coord$lat_2 <- as.numeric(gsub("[^0-9.-]", "", aqui_coord$lat))
aqui_coord$long_2 <- as.numeric(aqui_coord$long)
aqui_coord$total <- aqui_coord$g + aqui_coord$m
aqui_coord$total_norm <- 3 + (aqui_coord$total - min(aqui_coord$total)) * 7/
  (max(aqui_coord$total) - min(aqui_coord$total))

mplot <- ggplot() +
  borders("world", colour = "gray70", fill = "gray90") +
  geom_scatterpie(
    data = aqui_coord,
    aes(x = long_2, y = lat_2, r= total_norm),
    cols = c("g","m"),
    pie_scale = 0.5
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  coord_fixed()

mplot
