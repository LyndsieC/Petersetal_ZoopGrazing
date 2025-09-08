# Collis, June 2024
# Plot field chl (for Fig 2)

# Libraries ----
library(tidyverse)

#### load field chl data #
field_chl <- read.csv("00_ZP_grazing_field_chl_WBLE_2021.csv") %>% 
  select(Time, insitu_chla_ugL) %>% 
  mutate(insitu_chla_ugC_L = insitu_chla_ugL * 30, # 30 is the carbon-chla conversion factor (determined from Strickland 1960 and Geider 1987)
         Event = case_when(Time == "EarlyAug" ~ "August",
                           Time == "EarlySep" ~ "September",
                           TRUE ~ "NA"),
         Type = "Chlorophyll-a",
         Phyto_grp2 = "Total phytoplankton",
         bms_ugCL = insitu_chla_ugC_L)

# plot field chl
png(file = "14_WBLE_FieldChla_Fige2b.png",
    width = 600, height = 800)

ggplot(field_chl, aes(Event, insitu_chla_ugC_L)) + 
  geom_col(just = 0.5) +
  scale_x_discrete(labels = c("August", "September")) +
  ylab(expression(paste("Phytoplankton Biomass (µg C  ", L^-1,")"))) +
  labs(tag = "b)") +
  theme_bw() +
  theme(text = element_text(family = "serif"),
        axis.title.y = element_text(size = 24, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 24, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 19),
        legend.position = "top",
        legend.title= element_blank(), 
        panel.background = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 20),
        strip.background.y = element_rect(fill = "transparent", color = "transparent", linewidth = 1),
        strip.background.x = element_rect(fill = "transparent")) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

dev.off()
