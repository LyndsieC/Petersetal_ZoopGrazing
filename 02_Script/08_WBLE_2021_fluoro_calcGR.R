# Collis, June 2023
# Calculate grazing rate, % SS consumed, % upot consumed for fluoroprobe phyto groups
# Doing only for mesozoops- Early Aug and Early Sept for the nano/pico grazing MS

# Libraries ----
library(tidyverse)
library(broom)

# load data ----
# Load zoop clearance rate and phyto growth rate data for extracted chl-a, calcuated in 05_WBLE_2021_chl_calcCRGR_boot_LMC
fluoro_CR_growrate <- read.csv("07_CR_growrate_est_2021_36873_ZPgrazing_fluoro_boot100.csv") %>% 
  select(-X) %>% 
  mutate(Time = as.factor(Time),
         Phyto_Grp = as.factor(Phyto_Grp))

#### load field fluoroprobe data to calculate grazing rates ####
fluoro_field <- read.csv("00_fluoroprobe_EarlyAugSept_field_WBLE2021.csv")
# subset fluoroprobe data from 200um sieved samples, Early August and Early September
fluoro_field_200 <- fluoro_field %>% 
  # mutate(Phyto_grp = replace(Phyto_grp, Phyto_grp == "Bluegreen", "Cyanobacteria")) %>% 
  filter(Treatment == "200um_SIEVED")

fluoro_field_200_comb <- fluoro_field_200 %>% 
  select(Time, Phyto_Grp, Avg_ugL) %>% 
  pivot_wider(names_from = Phyto_Grp, values_from = Avg_ugL) %>% 
  mutate(CryptoDiatoms = Cryptophyta + Diatoms) %>% 
  select(-c(Cryptophyta, Diatoms)) %>% 
  pivot_longer(cols = Bluegreen:CryptoDiatoms,
               names_to = "Phyto_Grp",
               values_to = "chl_ugL") 

# get C conversion factors
C_conv <- read.csv("00_fluoroprobe_Cconversions.csv")

# Load zooplankton density data
# ZP_dens_field_num_L is field zoop density from LEPAS-- density estimate contains all organisms >200um
# ZP_dens_exp_num_L is the average zoop density in the experiment- calculated from the average of the zoop density in each of the zoop treatments
# ZP_field <- read.csv("Data_for_R/WBLE_2021_ZPgrazingexps_ZPdens.csv") %>% 
#   mutate(Time = as.factor(Time))

# Load zooplankton biomass data-- from field site and avg bms during exps
Zp_bms <- read.csv("00_LEPAS_WBLE_2021_36873_FieldExp_BMS_ugCL.csv") %>% 
  select(-X)

# join zoop and fluoro field data for GR calcs
fluoro_CR_growj <- fluoro_CR_growrate %>% 
  left_join(fluoro_field_200_comb, by = c("Time", "Phyto_Grp")) %>% 
  left_join(C_conv, by = "Phyto_Grp") %>% 
  left_join(Zp_bms, by = "Time")

# convert chl concentration to carbon units-- using conversion factor 60 C : 1 chl-a
fluoro_CR_grow_Cunits <- fluoro_CR_growj %>% 
  mutate(chl_ugC_L = chl_ugL * C_conversion) # multiply by the carbon-chla conversion factor

fluoro_CR_GR <- fluoro_CR_grow_Cunits %>% 
  filter(Phyto_Grp != "Yellow Subst.") %>% 
  mutate(GR_ugC_L_d_med = CR_L_ugC_d_med * chl_ugC_L * Zp_bms_ugC_L,
         GR_ugC_L_d_10 = CR_L_ugC_d_10 * chl_ugC_L * Zp_bms_ugC_L,
         GR_ugC_L_d_90 = CR_L_ugC_d_90 * chl_ugC_L * Zp_bms_ugC_L,
         per_SS_grazed_med = (GR_ugC_L_d_med / chl_ugC_L) * 100,
         per_SS_grazed_10 = (GR_ugC_L_d_10 / chl_ugC_L) * 100,
         per_SS_grazed_90 = (GR_ugC_L_d_90 / chl_ugC_L) * 100,
         m_med = CR_L_ugC_d_med * Zp_bms_exp_ugC_L, # calculate grazing mortality rate- ZP dens is the mean ZP dens in the grazing experiment treatments
         m_10 = CR_L_ugC_d_10 * Zp_bms_exp_ugC_L,
         m_90 = CR_L_ugC_d_90 * Zp_bms_exp_ugC_L,
         m_fld_med = CR_L_ugC_d_med * Zp_bms_ugC_L, # calculate field grazing mortality rate-- ZP bms is from the field site
         m_fld_10 = CR_L_ugC_d_10 * Zp_bms_ugC_L,
         m_fld_90 = CR_L_ugC_d_90 * Zp_bms_ugC_L,
         per_upot_med = (m_med / Int_med) * 100,
         per_upot_10 = (m_10 / Int_med) * 100, # using median Intercept (u) for upper and lower bounds b/c that's the average food growth rate over the experiment
         per_upot_90 = (m_90 / Int_med) * 100,
         per_upot_fld_med = (m_fld_med / Int_med) * 100,
         per_upot_fld_10 = (m_fld_10 / Int_med) * 100, # using median Intercept (u) for upper and lower bounds b/c that's the average food growth rate over the experiment
         per_upot_fld_90 = (m_fld_90 / Int_med) * 100)

# export raw values to csv
write.csv(fluoro_CR_GR, "08_2021_36873_CRGR_fluoro_raw.csv")

# When the lower 10th percentile for the ind CR is <0, this indicates there wasn't significant grazing, so make the GR and GR zero
fluoro_CR_GR2 <- fluoro_CR_GR %>% 
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
         per_upot_90 = if_else(CR_mL_ugC_d_10 > 0, per_upot_90, 0),
         per_upot_fld_med = if_else(CR_mL_ugC_d_10 > 0, per_upot_fld_med, 0),
         per_upot_fld_10 = if_else(CR_mL_ugC_d_10 > 0, per_upot_fld_10, 0),
         per_upot_fld_90 = if_else(CR_mL_ugC_d_10 > 0, per_upot_fld_90, 0))

# export "cleaned-up" data (w/ zeros) to csv
write.csv(fluoro_CR_GR2, "08_2021_36873_CRGR_fluoro.csv")
