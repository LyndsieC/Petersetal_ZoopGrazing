# Collis/Peters, March 2022
# Calculate biovolume of nano- and picoplankton from field samples collected during 2021 field season at LEPAS site 36-873
# Nano- and picoplankton counts were conducted by Dan Peters via epifluorescence microscopy

# load libraries
library(tidyverse)
library(viridis)
library(hrbrthemes)

# load nano and picoplankton counts (from epifluorescence microscopy)
IC_dat <- read.csv("NanoPico_InitialCollection_Count.csv")

# manipulate df
IC_dat2 <- IC_dat %>% 
  separate(Name,c("Treatment","Trial","InitialFinal","MesoMicro","Time","VolumeFiltered","DAPIconc","FilterSize", # create new columns from the original "Name" column,
                  "FilterType","ImageNum", "nd2"), "_") %>%  # each column name specified here
  mutate(Size = case_when(FilterSize == "0.2f" ~ "Pico", # create new column to differentiate Nano vs. Pico counts
                          FilterSize == "0.8f" ~ "Nano", # differentiate nano vs. pico based on filter size
                          TRUE ~ "NA"))

# Calculate biovolume of each individual particle
IC_dat2a <- IC_dat2 %>% 
  mutate(bv_um3 = ((Length * 3.14 * (Width / 2)^2) + ((4/3) * 3.14 * (Width/2)^3) - (Width * 3.14 * (Width/2)^2)))

# Clean up data and get sum of biovolume for each image
IC_dat3 <- IC_dat2a %>% 
  group_by(Time, VolumeFiltered, FilterSize, FilterType, ImageNum) %>% # summarize the data by the categories contained within the specified columns
  summarise(Count = n(), bv_um3 = sum(bv_um3)) # summarize by the number of cells counted for a given image in a given nano/pico sample

# remove "mL" from the volume filtered column
# Have to assign volume as "character" for now to get it to work, I make it numeric in the next biovolume equation
IC_dat3$VolumeFiltered <- gsub("mL","", as.character(IC_dat3$VolumeFiltered))

# get average biovolume and count across all the images
IC_dat4 <- IC_dat3 %>% 
  group_by(Time, VolumeFiltered, FilterSize, FilterType) %>% 
  summarise(avg_count = mean(Count),
            avg_bv_um3 = mean(bv_um3))

# convert average image BV to average filter BV (you can think of this as the "total BV on the filter", estimated from the average BV of the images)
# This is calculated as average image BV * (Filter Area (Af) / Image Area (Af))
# Filter area is calculated as r^2*pi, where r^2 is the radius of the filter
# Filter diameter is 25mm, so radius is 12.5mm or 12500um
Af <- (12500)^2 * pi
# Image area is calculated as image length x width
# Image length is 82.56um, Image width is 66.05um
Ai <- 66.05 * 82.56
# calculate the estimated biovolume on the filter
IC_dat5 <- IC_dat4 %>% 
  mutate(filterBV_um3 = avg_bv_um3 * (Af / Ai))

# convert to average BV per mL for each sample/filter type
IC_dat6 <- IC_dat5 %>% 
  mutate(avgBV_um3_mL = filterBV_um3 * (1 / as.numeric(VolumeFiltered)),
         Time = fct_relevel(Time, 
                             "LateApril", "LateMay", "LateJune", "MidJuly","LateJuly", "EarlyAug", 
                             "LateAug", "EarlySep", "LateSep"))

ggplot(IC_dat6, aes(y = avgBV_um3_mL, color = FilterSize, x = Time))+
  geom_point()+
  # geom_smooth(method = "lm",se = FALSE)+
  facet_grid(FilterSize ~ FilterType, scales = "free")

# pivot wider by DAPI and Chl to calculate hetertroph BV
IC_dat7 <- IC_dat6 %>% 
  pivot_wider(names_from = FilterType, 
              values_from = c(avg_count, avg_bv_um3, filterBV_um3, avgBV_um3_mL))

# calculate heterotroph BV
IC_dat8 <- IC_dat7 %>% 
  mutate(avgBV_um3_mL_Hetero = avgBV_um3_mL_Dapi - avgBV_um3_mL_Chl)

# convert data back to long format
IC_dat9 <- IC_dat8 %>% 
  select(Time, FilterSize, avgBV_um3_mL_Chl, avgBV_um3_mL_Dapi, avgBV_um3_mL_Hetero) %>% 
  pivot_longer(cols = avgBV_um3_mL_Chl:avgBV_um3_mL_Hetero, names_to = "FilterType", values_to = "IC_BV_um3_mL") %>% 
  separate(FilterType, into = c("a", "b", "c", "FilterType"), sep = "_") %>% 
  select(-c(a, b, c)) %>% 
  mutate(IC_BV_um3_L = IC_BV_um3_mL * 1000) # convert units to L

# export to csv
write.csv(IC_dat9, "03_WBLE_2021_InitialCollection_NanoPico_BV.csv")


## Create plot for MS-- Early August and Early September
# make dataframe for plotting
plot_df <- IC_dat9 %>% 
  filter(Time == "EarlyAug" | Time == "EarlySep") %>% 
  filter(FilterType != "Dapi") %>% 
  mutate(Date = case_when(Time == "EarlyAug" ~ "August",
                          Time == "EarlySep" ~ "September"),
         Size = case_when(FilterSize == "0.2filter" ~ "Picoplankton",
                          FilterSize == "0.8filter" ~ "Nanoplankton"),
         Troph = case_when(FilterType == "Chl" ~ "Autotroph",
                           FilterType == "Hetero" ~ "Heterotroph"),
         bms_ugC_L = (IC_BV_um3_L * 0.36) / 1000000, # Verity et al. 1992 C-bv conversion: 0.36 pg/um^3
         Size = factor(Size, levels = c("Picoplankton", "Nanoplankton")))
         
# plot nano/pico initial collection for MS-- Early August and Early Sept 
png(file = "03_IC_NanoPico_biomass_Fig3.png",
    width = 510, height = 500)
    
plot_df %>% 
  filter(Troph != "Total") %>% 
  ggplot(aes(y = bms_ugC_L, x = Date, fill = Troph)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(.~Size) +
  labs(y = expression(paste("Biomass  (µg C  ",  L^{-1},")"))) +
  # scale_fill_viridis(discrete = T) +
  # scale_fill_manual(values = c("darkorchid4", "darkgreen")) +
  scale_fill_manual(values = c("darkgreen", "#F0E442")) +
  theme_bw() +
  theme(text = element_text(family = "serif"),
    axis.text.x = element_text(size = 24, color = "black"),
    axis.text.y = element_text(size = 21, color = "black"),
    axis.title.y = element_text(size = 26, color = "black"),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="top",
    legend.title = element_blank(),
    legend.text=element_text(size = 20),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 23),
    panel.spacing = unit(0,'lines'),
    plot.tag = element_text(size = 17),
    plot.tag.position = c(0.13, 0.69))

dev.off()
