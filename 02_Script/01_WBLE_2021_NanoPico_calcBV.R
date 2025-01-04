# Peters/Collis, December 2021
# Updated March 2022
# Calculate biovolume of nano- and picoplankton from zooplankton grazing experiments (Aug and Sept 2021)
# Nano- and picoplankton counts were conducted by Dan Peters via epifluorescence microscopy

# load libraries
library(tidyverse)
library(cowplot)
library(stringr)

## Load data ####
# read csv files into R-- individual nano/pico measurements from epifluorescence microscopy
# Length and Width are in micrometer (um) units
nanoEAug <- read.csv("NanoEAugCount.csv")
nanoESep <- read.csv("NanoESepCount.csv")
picoEAug <- read.csv("PicoEAugCount.csv")
picoESep <- read.csv("PicoESepCount.csv", stringsAsFactors=FALSE, fileEncoding="latin1")
# combine all files
alldat <- rbind(nanoEAug, nanoESep, picoEAug, picoESep)

# load sample info
sampinfo_forbv <- read.csv("00_NanoPico_WBLE_2021_SampleInfo.csv") %>% 
  select(-c(InitialTime, FinalTime, Density_numL)) %>%  # remove these columns to help with bv estimate and df manipulation
  mutate(Treatment = as.numeric(Treatment))

# check data structure
str(alldat)

# create new columns from the original "Name" column
alldat2 <- alldat %>% 
  separate(Name, c("Treatment", "Trial", "InitialFinal", "MesoMicro", "Time", "VolumeFiltered", "DAPIconc", "FilterSize", 
                  "FilterType", "ImageNum", "nd2"),"_") %>% # each column name specified here
  mutate(Size = case_when(FilterSize == "0.2filter" ~ "Pico", # create new column to differentiate Nano vs. Pico counts
                          FilterSize == "0.8filter" ~ "Nano", # differentiate nano vs. pico based on filter size
                          TRUE ~ "NA"))

## This section is a "check" on the data ####
# get total number of organisms counted in each image
alldat_ct <- alldat2 %>% 
  group_by(Time, Size, Treatment, Trial, InitialFinal, FilterType, ImageNum) %>% 
  summarise(Count = n()) # summarize by the number of cells counted for a given image in a given nano/pico sample

# plot organism counts
ggplot(alldat_ct %>% filter(Size == "Nano"),aes(y = Count, x = Treatment, color = InitialFinal, shape = Trial)) +
  geom_point() +
  facet_grid(Time ~ FilterType)

ggplot(alldat_ct %>% filter(Size == "Pico"),aes(y = Count, x = Treatment, color = InitialFinal, shape = Trial)) +
  geom_point() +
  facet_grid(Time ~ FilterType)

# get average count per sample
alldat_avgct <- alldat_ct %>% 
  group_by(Time, Size, Treatment, Trial, InitialFinal, FilterType) %>% 
  summarise(avg_ct = mean(Count))

ggplot(alldat_avgct %>% filter(Size == "Nano"), aes(y = avg_ct, x = Treatment, color = InitialFinal, shape = Trial))+
  geom_point() +
  facet_grid(Time ~ FilterType)

ggplot(alldat_avgct %>% filter(Size == "Pico"), aes(y = avg_ct, x = Treatment, color = InitialFinal, shape = Trial))+
  geom_point() +
  facet_grid(Time ~ FilterType)

# spread data by initial/final to calculate change in nano/pico over the 24 hour experiment
alldats <- alldat_avgct %>% 
  spread(InitialFinal, avg_ct) %>% 
  mutate(f_i = Final / Initial)

# plot f/i for the experiments
ggplot(alldats %>% filter(Size == "Nano"), aes(y=f_i, x = Treatment, shape = Trial))+
  geom_point()+
  facet_grid(Time~FilterType)

ggplot(alldats %>% filter(Size == "Pico"), aes(y=f_i, x = Treatment, shape = Trial))+
  geom_point()+
  facet_grid(Time~FilterType)

## Biovolume calculations ####
# calculate biovolume of each organism (um^3), as per Smayda 1978
alldat2a <- alldat2 %>% 
  select("Time", "Size", "Treatment", "Trial", "InitialFinal", "ImageNum", "FilterType", "Length", "Width", "Area", "VolumeFiltered") %>% 
  mutate(bv_um3 = ((Length * 3.14 * (Width / 2)^2) + ((4/3) * 3.14 * (Width/2)^3) - (Width * 3.14 * (Width/2)^2)))

# get total count and biovolume for each image
alldat3 <- alldat2a %>% 
  group_by(Time, Size, Treatment, Trial, InitialFinal, ImageNum, FilterType) %>% # summarize the data by the categories contained within the specified columns
  summarise(Count = n(), # sum the number of cells counted for a given image in a given nano/pico sample
            bv_um3 = sum(bv_um3, na.rm = TRUE)) # sum the biovolume for a given image in a given nano/pico sample

ggplot(alldat3 %>% filter(Size == "Nano"), aes(y = bv_um3, x = Treatment, color = InitialFinal, shape = Trial))+
  geom_point()+
  facet_grid(Time ~ FilterType)

ggplot(alldat3 %>% filter(Size == "Pico"), aes(y = bv_um3, x = Treatment, color = InitialFinal, shape = Trial))+
  geom_point()+
  facet_grid(Time ~ FilterType)

ggplot(alldat3 %>% filter(Size == "Nano"), aes(y = bv_um3, x = Treatment, color = InitialFinal)) +
  geom_boxplot() +
  facet_grid(Time~FilterType)

# "spread" data so chl and dapi counts are separate columns
alldat4 <- alldat3 %>% 
 select(-Count) %>% 
  spread(FilterType, bv_um3)

# Calculate heterotroph biovolume
# We're doing this by subtracting DAPI biovolume from chl biovolume-- this is because you can't use the individual cells to figure out which are heterotrophs-- you have to make the calculation at the "population" level for each image
# Some heterotroph values are negative because estimated chl BV was higher than DAPI BV-- this occurred because some images had high fluorescence for chl, so the chl BV was higher than the DAPI estimated BV
alldat5 <- alldat4 %>% 
  mutate(Hetero = Dapi - Chl) # Heterotroph biovolume = Total biovolume (DAPI) - Autotrophic biovolume (Chl)

# get average biovolume of each filtered sample (average of all images)
alldat6 <- alldat5 %>% 
  group_by(Time, Size, Treatment, Trial, InitialFinal) %>% 
  summarize(Chl = mean(Chl, na.rm = TRUE),
            Dapi = mean(Dapi, na.rm = TRUE),
            Hetero = mean(Hetero, na.rm = TRUE))

# convert df to long format
alldat7 <- alldat6 %>% 
  pivot_longer(
    cols = Chl:Hetero,
    names_to = "FilterType",
    values_to = "bv_um3") %>% 
  mutate(Treatment = as.numeric(Treatment))

## Convert biovolume from (um^3) -> (um^3) / mL ####
# to do this, we're taking the average biovolume on the filter (obtained by averaging all image biovolumes), and accounting for the volume of 
# water filtered onto the nucleopore filter (1 mL for pico, 10 mL for nano), as well as the total area of the filter vs. the area of each image

# join sample info df to get volume filtered
alldatj <- alldat7 %>% 
  left_join(sampinfo_forbv, by = c("Time", "Size", "Treatment", "Trial"))

# Filter calculations to convert biovolume from (um^3) -> (um^3) / mL
# Microscope photo image area (L x W):  82.56 um * 66.05 um = 5453.088 um^2
# Convert um^2 to mm^2: 5453.088 / 1000000 = 0.005453088 mm^2
# Filter area (25mm filter-- pi x r^2): 3.1415 x 12.5^2 mm = 490.8592 mm^2
# Unitless conversion factor to account for filter area imaged: 0.005453088 / 490.8592 = 0.00001110928

# The mutate here is to adjust for area of filter and Volume filtered
# Converts units from (um^3) -> (um^3)/mL
alldat8 <- alldatj %>% 
  mutate(bv_um3_mL = (bv_um3 / 0.00001110928) / Volume_mL) # unitless conversion factor calculated above

alldat9 <- alldat8 %>% 
  select(-c(bv_um3, Volume_mL)) %>% 
  pivot_wider(
    names_from = InitialFinal,
    values_from = bv_um3_mL
  )

# export to csv
write.csv(alldat9, "01_WBLE_2021_NanoPico_biovolume_um3_mL.csv")
