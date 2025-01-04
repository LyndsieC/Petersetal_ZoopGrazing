# Collis, June 2023
# Updated June 2024 for biomass-corrected data
# Combine nano/pico, chl, and fluoroprobe CR/GR rates into one dataframe

# Libraries ----
library(tidyverse)
library(broom)

# load data ----
# nano/picoplankton
nanopico_dat <- read.csv("04_2021_36873_CRGR_bms_NanoPico.csv") %>% 
  select(-X) %>% 
  mutate(Date = case_when(Time == "EarlyAug" ~ "August",
                          Time == "EarlySep" ~ "September",
                          TRUE ~ "blah"),
         Facet = case_when(Size == "Nano" ~ "Nanoplankton",
                           Size == "Pico" ~ "Picoplankton",
                           TRUE ~ "blah"),
         FoodType = case_when(Filter == "Dapi" & Size == "Nano" ~ "Total Nanoplankton",
                              Filter == "Chl" & Size == "Nano" ~ "Autotrophic Nanoplankton",
                              Filter == "Hetero" & Size == "Nano" ~ "Heterotrophic Nanoplankton",
                              Filter == "Dapi" & Size == "Pico" ~ "Total Picoplankton",
                              Filter == "Chl" & Size == "Pico" ~ "Autotrophic Picoplankton",
                              Filter == "Hetero" & Size == "Pico" ~ "Heterotrophic Picoplankton",
                              TRUE ~ "blah"),
         FoodAbund_ugC_L = IC_BV_ugC_L)

# extracted chla
chl_dat <- read.csv("06_2021_36873_CRGR_chl.csv") %>% 
  select(-X) %>% 
  mutate(Date = case_when(Time == "EarlyAug" ~ "August",
                          Time == "EarlySep" ~ "September",
                          TRUE ~ "blah"),
         Facet = "Chlorophyll-a",
         FoodType = "Extracted Chl-a",
         FoodAbund_ugC_L = insitu_chla_ugC_L)

# fluoroprobe data
fluoro_dat <- read.csv("08_2021_36873_CRGR_fluoro.csv") %>% 
  select(-X) %>% 
  mutate(Date = case_when(Time == "EarlyAug" ~ "August",
                          Time == "EarlySep" ~ "September",
                          TRUE ~ "blah"),
         Facet = "Fluoroprobe",
         FoodType = case_when(Phyto_Grp == "CryptoDiatoms" ~ "Cryptophyta and Bacillariophyta",
                              Phyto_Grp == "Bluegreen" ~ "Cyanobacteria",
                              Phyto_Grp == "Green Algae" ~ "Chlorophyta",
                              TRUE ~ "blah"),
         FoodAbund_ugC_L = chl_ugC_L)

all_dat <- bind_rows(nanopico_dat, chl_dat, fluoro_dat)

all_dat2 <- all_dat %>%
  mutate(FoodType = fct_relevel(FoodType, 
                                "Total Nanoplankton", 
                                "Autotrophic Nanoplankton", 
                                "Heterotrophic Nanoplankton",
                                "Total Picoplankton",
                                "Autotrophic Picoplankton",
                                "Heterotrophic Picoplankton",
                                "Cryptophyta and Bacillariophyta",
                                "Cyanobacteria",
                                "Chlorophyta",
                                "Extracted Chl-a"),
         Date = fct_relevel(Date, 
                            "August", "September"))

# get max CR, GR, and food abundance for each date
max_dat <- all_dat2 %>% 
  group_by(Time) %>% 
  summarize(GR_max = max(GR_ugC_L_d_med),
            CR_max = max(CR_mL_ugC_d_med),
            Food_max = max(FoodAbund_ugC_L))

# join max values to data
all_dat_max <- all_dat2 %>% 
  left_join(max_dat, by = "Time")

# calculate relative clearance rate, grazing rate, and relative food abundance
all_dat3 <- all_dat_max %>% 
  mutate(CR_Ratio_Med = CR_mL_ugC_d_med / CR_max,
         GR_Ratio_Med = GR_ugC_L_d_med / GR_max,
         Foodi_FoodMax = FoodAbund_ugC_L / Food_max)

# reorder factors
all_dat4 <- all_dat3 %>%
  mutate(FoodType = fct_relevel(FoodType, 
                                "Extracted Chl-a", "Cryptophyta and Bacillariophyta",
                                "Chlorophyta", "Cyanobacteria",
                                "Total Nanoplankton","Autotrophic Nanoplankton", 
                                "Heterotrophic Nanoplankton","Total Picoplankton", 
                                "Autotrophic Picoplankton", "Heterotrophic Picoplankton"),
         Date = fct_relevel(Date, "August", "September"))

all_dat_fin <- all_dat4 %>% 
  select(Time, Date, Int_med:CR_L_ugC_d_90, Zp_bms_ugC_L, Zp_bms_exp_ugC_L, GR_ugC_L_d_med:per_upot_90, Facet, FoodType, FoodAbund_ugC_L, GR_max:Foodi_FoodMax)

# export data to csv
write.csv(all_dat_fin, "09_2021_36873_NanoPicoMS_CRGR_final.csv")
