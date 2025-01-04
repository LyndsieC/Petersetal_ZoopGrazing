# load libraries
library(tidyverse)

# load files
# LEPAS phyto-- bms less than 200um calculated in 01_LEPAS_phyto_lengths_NanoPicoMS_031724.R
lepas <- read.csv("15a_LEPAS_phyto_bms_less200um_2021_36873.csv")
# fluoroprobe-- bms less than 200um (lake water was sieved through a 200um sieve for the grazing experiments)
fluoroprobe <- read.csv("00_fluoro_phyto_bms_ugCL_less200um_2021_36873.csv") 

# clean up data
lepas2 <- lepas %>% 
  select(-X) %>% 
  select(c(Sample_date, Phylum, bms_ugC_less200_tot)) %>% 
  rename(Date = Sample_date,
         Phyto_grp = Phylum,
         bms_ugCL = bms_ugC_less200_tot) %>% 
  mutate(Type = "Microscope")

fluoroprobe2 <- fluoroprobe %>% 
  select(-X) %>% 
  select(c(SampleDate, Phyto_grp2, C_ugL)) %>% 
  rename(Date = SampleDate,
         Phyto_grp = Phyto_grp2,
         bms_ugCL = C_ugL) %>% 
  mutate(Type = "FluoroProbe")

comb_dat <- rbind(lepas2, fluoroprobe2) %>% 
  mutate(Event = case_when(Date == "2021-08-09" | Date == "8/9/2021" ~ "August",
                           Date == "2021-09-10" | Date == "9/10/2021" ~ "September",
                           TRUE ~ "NA")) %>% 
  mutate(Phyto_grp2 = factor(Phyto_grp, 
                             levels = c("Cryptophyta and Bacillariophyta", "Cryptophyta", "Bacillariophyta", "Chlorophyta", "Cyanobacteria")))

# plot together
png(file = "15b_LEPASflouro_comb_less200um_C_ugL_2021_36873_Fig2a.png",
    width = 800, height = 700)

ggplot(comb_dat, aes(y = bms_ugCL, x = Event, fill = Phyto_grp2)) +
  geom_bar(stat = "identity") +
  facet_wrap(vars(Type), scales = "free_y") + # this allows the y-axis to be scaled differently for fluoroprobe vs. microscope
  scale_fill_manual(values = c("Cyanobacteria" = "deepskyblue3",
                               "Chlorophyta" = "chartreuse3",
                               "Cryptophyta" = "darkgoldenrod1",
                               "Bacillariophyta" = "darkgoldenrod3",
                               "Cryptophyta and Bacillariophyta" = "darkred")) +
  scale_x_discrete(labels = c("August", "September")) +
  ylab(expression(paste("Phytoplankton Biomass (µg C  ", L^-1,")"))) +
  labs(tag = "a)") +
  theme_bw() +
  theme(text = element_text(family = "serif"),
        axis.title.y = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 23),
        legend.position = "top",
        legend.title= element_blank(), 
        panel.background = element_rect(fill = "transparent", linewidth = 1),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 25),
        strip.background.y = element_rect(fill = "transparent", color = "transparent", linewidth = 1),
        strip.background.x = element_rect(fill = "transparent"),
        plot.tag = element_text(size = 16),) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

dev.off()
