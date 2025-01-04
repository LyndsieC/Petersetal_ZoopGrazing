# load libraries
library(tidyverse)
library(lubridate)
library(scales)
library(devtools)
library(ggh4x)

# load data
fluoro <- read.csv("InitialCollection_2021_ZPexps_rawdata4.csv")

fluorow <- fluoro %>%
  mutate(Date_Time = as.POSIXct(Date_Time, format = "%m/%d/%Y %H:%M")) %>%
  mutate(Bac_Cryp = Bacillariophyta + Cryptophyta,
         Total = Bacillariophyta + Cryptophyta + Cyanobacteria + Cryptophyta)

fluorol <- fluorow %>%
  pivot_longer(Chlorophyta:Total, names_to = 'Taxa', values_to = 'Conc_ugL')

#### Create Appendix Figure 3-- all phyto groups #####################
fluorol2 <- fluorol %>% 
  filter(Taxa != "YellowSubst")

#### High HABs (Crypots and Diatoms interference)-- Late Sept, rep 2
late_sept <- fluorol2 %>% 
  filter(Type == "SIEVED_200um") %>% 
  filter(Event == "Late_September" & Rep == "2") %>% 
  group_by(Taxa) %>% 
  mutate(Sec = seq(1, 77, by = 2),  # create column of time elapsed for each measurement (fluoroprobe readings taken every 2 seconds
         Taxa = factor(Taxa, levels = c("Chlorophyta", "Cyanobacteria", "Bacillariophyta", "Cryptophyta", "Bac_Cryp", "Total")), # reorder factors
         Time = "High cHAB")

#### Low HABs (Crypots and Diatoms no interference)-- Late June, rep 2
late_june <- fluorol2 %>% 
  filter(Type == "SIEVED_200um") %>% 
  filter(Event == "Late_June" & Rep == "2") %>% 
  group_by(Taxa) %>% 
  mutate(Sec = seq(1, 77, by = 2),  # create column of time elapsed for each measurement (fluoroprobe readings taken every 2 seconds
         Taxa = factor(Taxa, levels = c("Chlorophyta", "Cyanobacteria", "Bacillariophyta", "Cryptophyta", "Bac_Cryp", "Total")), # reorder factors
         Time = "Low/No cHAB")

# Plot both dates together to compare
fluoro_comp <- rbind(late_sept, late_june)

png(file = "13_WBLE_Fluoroprobe_SuppFige3.png",
    width = 1000, height = 800)

ggplot(fluoro_comp %>% filter(Sec <= 61), aes(y = Conc_ugL, x = Sec, color = Taxa)) +
  geom_line(size = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(labels = c("Chlorophyta", "Cyanobacteria", "Bacillariophyta", "Cryptophyta", "Cryptophyta +\nBacillariophyta", "Total"),
                     values = c("Cyanobacteria" = "deepskyblue3",
                                "Cryptophyta" = "darkgoldenrod1",
                                "Chlorophyta" = "chartreuse3",
                                "Bacillariophyta" = "darkred",
                                "Bac_Cryp" = "darkgrey",
                                "Total" = "black")) +
  labs( x = "Elapsed time (sec)", 
        y = expression(paste("Chlorophyll  (µg  ",  L^{-1},")"))) +
  ggh4x::facet_grid2(. ~ Time, scales = "free_y", independent = "y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=16, color = "black"),
    axis.text.y = element_text(size=16, color = "black"),
    axis.title.y = element_text(size=18),
    axis.title.x = element_text(size=18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x=unit(0.5, "lines"),
    legend.position="top",
    legend.title = element_blank(),
    legend.text=element_text(size=12),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 16),
    panel.spacing = unit(0,'lines'),
    plot.tag = element_text(size = 16),
    plot.tag.position = c(0.13, 0.69))

dev.off()
