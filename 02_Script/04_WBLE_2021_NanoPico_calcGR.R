# Collis, March 2023
# Updated June 2024 to change units to zooplankton biomass
# Calculate grazing rates of nano and picoplankton, based on clearance rates obtained in script 03
# Doing only for mesozoops- Early Aug and Early Sept for the nano/pico grazing MS

# load libraries
library(tidyverse)
library(cowplot)
library(stringr)
library(broom)

# Load zoop clearance rate and phyto growth rate data, calcuated in 02_WBLE_2021_NanoPico_calcCRGR_boot_LMC
CR_growrate <- read.csv("02_CR_growrate_bms_2021_36873_ZPgrazing_boot100.csv") %>% 
  select(-X) %>% 
  mutate(Time = as.factor(Time),
         Size = as.factor(Size),
         Filter = as.factor(Filter))

# Load "Initial Collection" (i.e., field collection) nano- and picoplankton biovolume estimates, calculated in 03_InitialCollection_NanoPico_calcBV_LMC
IC_bv <- read.csv("03_WBLE_2021_InitialCollection_NanoPico_BV.csv") %>% 
  mutate(Size = case_when(FilterSize == "0.2filter" ~ "Pico",
                            FilterSize == "0.8filter" ~ "Nano",
                            TRUE ~ "NA")) %>% 
  rename(Filter = FilterType) %>% 
  mutate(Time = as.factor(Time),
         Size = as.factor(Size),
         Filter = as.factor(Filter)) %>% 
  select(Time, Size, Filter, IC_BV_um3_mL, IC_BV_um3_L)

# Load zooplankton biomass data-- from field site and avg bms during exps
Zp_bms <- read.csv("00_LEPAS_WBLE_2021_36873_FieldExp_BMS_ugCL.csv") %>% 
  select(-X)

#### Calculate grazing rate ####
# Join CR/growth rate df w/ Initial Collection nano/pioc biovolume and mesozoop density
CR_growratej <- CR_growrate %>% 
  left_join(IC_bv, by = c("Time", "Size", "Filter")) %>% 
  left_join(Zp_bms, by = "Time")

# convert nano/pico biovolume to carbon units-- using conversion factor 0.36pgC/um3 (Verity et al. 1992)
CR_growrate_Cunits <- CR_growratej %>% 
  mutate(IC_BV_ugC_L = (IC_BV_um3_L * 0.36) / 1000000) # dividing by 1000000 converts pg to ug, 0.36 is the carbon conversion factor

CR_GR <- CR_growrate_Cunits %>% 
  mutate(CR_for_perPP_med = ifelse(CR_L_ugC_d_10 < 0, 0, CR_L_ugC_d_med), # to calculate % prim prod grazed, I need to force neg CR values to zero so the m term isn't negative
         CR_for_perPP_10 = ifelse(CR_L_ugC_d_10 < 0, 0, CR_L_ugC_d_10),
         CR_for_perPP_90 = ifelse(CR_L_ugC_d_10 < 0, 0, CR_L_ugC_d_90)) %>% 
  mutate(GR_ugC_L_d_med = CR_L_ugC_d_med * IC_BV_ugC_L * Zp_bms_ugC_L,
         GR_ugC_L_d_10 = CR_L_ugC_d_10 * IC_BV_ugC_L * Zp_bms_ugC_L,
         GR_ugC_L_d_90 = CR_L_ugC_d_90 * IC_BV_ugC_L * Zp_bms_ugC_L,
         per_SS_grazed_med = (GR_ugC_L_d_med / IC_BV_ugC_L) * 100,
         per_SS_grazed_10 = (GR_ugC_L_d_10 / IC_BV_ugC_L) * 100,
         per_SS_grazed_90 = (GR_ugC_L_d_90 / IC_BV_ugC_L) * 100,
         m_med = CR_L_ugC_d_med * Zp_bms_exp_ugC_L, # calculate grazing mortality rate- ZP bms is the mean ZP bms in the grazing experiment treatments
         m_10 = CR_L_ugC_d_10 * Zp_bms_exp_ugC_L,
         m_90 = CR_L_ugC_d_90 * Zp_bms_exp_ugC_L,
         per_upot_med = (m_med / Int_med) * 100,
         per_upot_10 = (m_10 / Int_med) * 100, # using median Intercept (u) for upper and lower bounds b/c that's the average food growth rate over the experiment
         per_upot_90 = (m_90 / Int_med) * 100)

# export raw values to csv
write.csv(CR_GR, "04_2021_36873_CRGR_bms_NanoPico_raw.csv")

# When the lower 10th percentile for the ind CR is <0, this indicates there wasn't significant grazing, so make the CR and GR zero
CR_GR2 <- CR_GR %>% 
  mutate(CR_mL_ugC_d_med = if_else(CR_mL_ugC_d_10 > 0, CR_mL_ugC_d_med, 0),
         CR_mL_ugC_d_10 = if_else(CR_mL_ugC_d_10 > 0, CR_mL_ugC_d_10, 0),
         CR_mL_ugC_d_90 = if_else(CR_mL_ugC_d_10 > 0, CR_mL_ugC_d_90, 0),
         GR_ugC_L_d_med = if_else(CR_mL_ugC_d_10 > 0, GR_ugC_L_d_med, 0),
         GR_ugC_L_d_10 = if_else(CR_mL_ugC_d_10 > 0, GR_ugC_L_d_10, 0),
         GR_ugC_L_d_90 = if_else(CR_mL_ugC_d_10 > 0, GR_ugC_L_d_90, 0),
         per_SS_grazed_med = if_else(CR_mL_ugC_d_10 > 0, per_SS_grazed_med, 0),
         per_SS_grazed_10 = if_else(CR_mL_ugC_d_10 > 0, per_SS_grazed_10, 0),
         per_SS_grazed_90 = if_else(CR_mL_ugC_d_10 > 0, per_SS_grazed_90, 0),
         per_upot_med = if_else(CR_mL_ugC_d_10 > 0, per_upot_med, 0),
         per_upot_10 = if_else(CR_mL_ugC_d_10 > 0, per_upot_10, 0),
         per_upot_90 = if_else(CR_mL_ugC_d_10 > 0, per_upot_90, 0))

# export "cleaned up" data (w/ zeros) to csv
write.csv(CR_GR2, "04_2021_36873_CRGR_bms_NanoPico.csv")
