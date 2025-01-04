## SET SYSTEM TIME 
Sys.setenv(TZ = "EST")

## LOAD LIBRARIES 
library(chron)
library(lubridate)
library(tidyverse)
library(gridExtra)  
library(ggthemes)
library(patchwork)

## Set theme for graphing
theme_set(theme_few(base_size = 16, base_family = "Helvetica"))
windowsFonts(Helvetica=windowsFont("TT Helvetica"))

# Load data
field12_18 <- read.csv("00_lake_erie_habs_field_sampling_results_2012_2018_v2.csv", check.names = F)
str(field12_18)
field19 <- read.csv("00_lake_erie_habs_field_sampling_results_2019.csv", check.names = F)
str(field19)

PAR_phyto_12_18 <- field12_18 %>% 
  select(c(Date:Sample_Depth_category, CTD_PAR_uE_m2_s, Extracted_Phycocyanin_ug_L, Extracted_Chla_ug_L)) %>% 
  mutate(Extracted_Phycocyanin_ug_L = as.numeric(Extracted_Phycocyanin_ug_L))
str(PAR_phyto_12_18)

PAR_phyto_19 <- field19 %>% 
  select(c(Date:Sample_Depth_category, CTD_PAR_uE_m2_s, Extracted_Phycocyanin_ug_L, Extracted_Chla_ug_L))
str(PAR_phyto_19)

field_PAR_phyto <- bind_rows(PAR_phyto_12_18, PAR_phyto_19) %>% 
  mutate(Date = as.POSIXct(Date, format = "%m/%d/%Y"),
         Month = month(Date),
         Year = year(Date),
         Day = yday(Date))

# get "water column" average (mixing in WBLE, so could expect phyto to experience an average light of the whole water column)
PAR_avg_WC <- field_PAR_phyto %>% 
  filter(Site == "WE2") %>% 
  filter(Sample_Depth_category != "Scum") %>% 
  summarize(PAR_avg = mean(CTD_PAR_uE_m2_s, na.rm = TRUE))

# subset surface samples
surface <- field_PAR_phyto %>% 
  filter(Sample_Depth_category == "Surface")

# subset WE2 site
WE2_surface <- surface %>% 
  filter(Site == "WE2")

# chla across season
# change label to month, but keep day-specific data plotted
p1 <- WE2_surface %>% 
  ggplot(aes(y = Extracted_Chla_ug_L, x = Day)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "loess", color = "black") + 
  annotate("text", x = 221, y = 51, label = "X", size = 8, color = "red") + # early august
  annotate("text", x = 253, y = 45, label = "X", size = 8, color = "red") + # early september date
  labs(tag = "a)") +
  ylab(expression(atop("Extracted chl-a",
                       paste("(µg ",L^-1,")")))) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=20, color = "black"),
    axis.title.y = element_text(size=20, color = "black"),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=12),
    legend.title = element_text(size = 13),
    legend.key.width= unit(1.5, 'cm'),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 14),
    panel.spacing = unit(0,'lines'),
    plot.tag = element_text(size = 16))

# phycocyanin across season
p2 <- WE2_surface %>%
  filter(Extracted_Phycocyanin_ug_L < 500) %>% 
  ggplot(aes(y = Extracted_Phycocyanin_ug_L, x = Day)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "loess", color = "black") +
  annotate("text", x = 221, y = 36, label = "X", size = 8, color = "red") + # early august date
  annotate("text", x = 253, y = 27, label = "X", size = 8, color = "red") + # early september date
  labs(tag = "b)") +
  ylab(expression(atop("Extracted phycocyanin",
                       paste("(µg ",L^-1,")")))) +
  xlab("Julian Day") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size=20, color = "black"),
        axis.title.y = element_text(size=20, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        legend.title= element_blank(), 
        legend.text=element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0,'lines'),
        plot.tag = element_text(size = 16))

# Combining chl and phycocyanin graphs
p12 <- p1 / p2

ggsave('12_WE2_HABscruise_chl_phyco_FigS2.png', p12,
       width = 8, height = 8)