# Collis, June 2023
# Calculate grazing rate, % SS consumed, % upot consumed for chl-a (extracted)
# Doing only for mesozoops- Early Aug and Early Sept for the nano/pico grazing MS

# Libraries ----
library(tidyverse)
library(broom)

# load data ----
# Load zoop clearance rate and phyto growth rate data for extracted chl-a, calcuated in 05_WBLE_2021_chl_calcCRGR_boot_LMC
chl_CR_growrate <- read.csv("05_CR_growrate_est_2021_36873_ZPgrazing_chl_bms_boot100.csv") %>% 
  select(-X) %>% 
  mutate(Time = as.factor(Time))

#### load field chl data to calculate grazing rates ####
field_chl <- read.csv("00_ZP_grazing_field_chl_WBLE_2021.csv") %>% 
  select(Time, insitu_chla_ugL) %>% 
  mutate(insitu_chla_ugC_L = insitu_chla_ugL * 30) # 30 is the carbon-chla conversion factor (determined from Strickland 1960 and Geider 1987)

# Load zooplankton biomass data-- from field site and avg bms during exps
Zp_bms <- read.csv("00_LEPAS_WBLE_2021_36873_FieldExp_BMS_ugCL.csv") %>% 
  select(-X)

#### Calculate grazing rate ####
# Join CR/growth rate df w/ Initial Collection nano/pioc biovolume and mesozoop density
chl_CR_growratej <- chl_CR_growrate %>% 
  left_join(field_chl, by = c("Time")) %>% 
  left_join(Zp_bms, by = "Time")

chl_CR_GR <- chl_CR_growratej %>% 
  mutate(GR_ugC_L_d_med = CR_L_ugC_d_med * insitu_chla_ugC_L * Zp_bms_ugC_L,
         GR_ugC_L_d_10 = CR_L_ugC_d_10 * insitu_chla_ugC_L * Zp_bms_ugC_L,
         GR_ugC_L_d_90 = CR_L_ugC_d_90 * insitu_chla_ugC_L * Zp_bms_ugC_L,
         per_SS_grazed_med = (GR_ugC_L_d_med / insitu_chla_ugC_L) * 100,
         per_SS_grazed_10 = (GR_ugC_L_d_10 / insitu_chla_ugC_L) * 100,
         per_SS_grazed_90 = (GR_ugC_L_d_90 / insitu_chla_ugC_L) * 100,
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
write.csv(chl_CR_GR, "06_2021_36873_CRGR_chl_raw.csv")

# When the lower 10th percentile for the ind CR is <0, this indicates there wasn't significant grazing, so make the GR and GR zero
chl_CR_GR2 <- chl_CR_GR %>% 
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

write.csv(chl_CR_GR2, "06_2021_36873_CRGR_chl.csv")
