# Collis/Peters/Hood, March 2022
# Updated May 2024, converted nano- and picoplankton to carbon units and expressed rates in terms of zooplankton biomass (carbon)
# Calculate bootstrapped estimates of indiv zoop clearance rates for nano- and picoplankton
# Doing only for mesozoops- Early Aug and Early Sept for the nano/pico grazing MS 

# load libraries
library(tidyverse)
library(cowplot)
library(stringr)
library(broom)

# load biovolume data, calculated in 01_WBLE_2021_NanoPico_calcBV.R
# units are um^3 / mL
alldat9 <- read.csv("01_WBLE_2021_NanoPico_biovolume_um3_mL.csv") %>% 
  select(-X)

# load sample info from experiments
sampinfo <- read.csv("00_NanoPico_WBLE_2021_SampleInfo.csv")

alldat9j <- alldat9 %>% 
  left_join(sampinfo, by = c("Time", "Size", "Treatment", "Trial")) 

ggplot(alldat9j, aes(y = Final, x = Initial)) +
  geom_point() +
  facet_grid(cols = vars(Time)) +
  geom_abline(intercept = 0)

# load field dens and bms-- will use to calculate avg ind dw and apply to zooplankton experiment treatments
zp_field_bms <- read.csv("00_LEPAS_WBLE_2021_36873_FieldExp_BMS_ugCL.csv") %>% 
  select(-X)

# calculate avg individual dry weight of zooplankton in each experiment, subset for August and September experiments only
ind_dw <- zp_field_bms %>% 
  filter(Time == "EarlyAug" | Time == "EarlySep") %>% 
  select(Time, Ind_wt_ugC)

#### Estimate clearance rates (CR)-- bootstrap CR estimates by randomly selecting 1 of the 4 Initial/Final bottle reps from each treatment ####
## Loop: calculate r (ln(Final/Initial)/Time) for a given Event/Size/ZP Treatment/Filter (chl, heterotroph, total-DAPI)-- bottle reps randomly generated
dat <- alldat9j

## create empty dataframe for bootstrapped data ----
datBoot <- as.data.frame(matrix(nrow = 1, ncol = 7))
names(datBoot) <- c("Time_i", "Filter_i", "Treatment_i", "Size_i", "FmI", "boot_i", "Rep_i")

# Time_i, Filter_i, Treatment_i, Size_i, FmI, boot_i, Rep_i

## Bootstrapping
# Bootstrap for 200 iterations
BootReps <- 100

# select event
for(e in 1:length(unique(dat$Time))){
  # Time_e <- unique(dat$Time)[1]
  Time_e <- unique(dat$Time)[e]
  dat_e <- dat %>% 
    filter(Time == Time_e)
  
  # select nano/pico size range
  for(s in 1:length(unique(dat$Size))){
    # Size_s <- unique(dat$Size)[1]
    Size_s <- unique(dat$Size)[s]
    dat_s <- dat_e %>% 
      filter(Size == Size_s)
    
    # select filter type
    for(f in 1:length(unique(dat$FilterType))){
      # Filter_f <-  unique(dat$FilterType)[1]
      Filter_f <- unique(dat$FilterType)[f]
      dat_f <- dat_s %>% 
        filter(FilterType == Filter_f)
      
      # select treatment
      for(t in 1:length(unique(dat_f$Treatment))){
        # Treatment_t <- unique(dat_f$Treatment)[1]
        Treatment_t <- unique(dat_f$Treatment)[t]
        dat_ft <- dat_f %>% 
          filter(Treatment == Treatment_t) %>% 
          mutate(FinalTime = as.POSIXct(FinalTime, format = "%m/%d/%Y %H:%M"),
                 InitialTime = as.POSIXct(InitialTime, format = "%m/%d/%Y %H:%M"))

        # within each filter-treatment combine all four finals with randomly selected initials
        for(i in 1:BootReps){
          # i = 1
          # log (Final/initial) where initials are randomly selected from four options within filter-treatment with replacement
          NumSeq_i <- sample(seq(1, 4, by = 1), size = 4, replace = 4)
        
          FmI <- log((dat_ft$Final/dat_ft[c(NumSeq_i),]$Initial)) / ((as.numeric(dat_ft$FinalTime) - as.numeric(dat_ft[c(NumSeq_i),]$InitialTime))) * 60 * 60 * 24
          
          # vector of event
          Time_i <- as.data.frame(rep(Time_e, 4))
          # rename event column
          names(Time_i) <- "Time_i"
          # vector of filter type
          Filter_i <- as.data.frame(rep(Filter_f, 4))
          # rename column in filter vector
          names(Filter_i) <- "Filter_i"
          # vector of rep numbers
          Rep_i <- dat_ft$Trial
          # vector of treatment values
          Treatment_i <- as.data.frame(rep(Treatment_t, 4))
          # rename treatment column
          names(Treatment_i) <- "Treatment_i"
          # vector of size values
          Size_i <- as.data.frame(rep(Size_s, 4))
          # rename size column
          names(Size_i) <- "Size_i"
          # boot strap number, ranges from 1 to BootReps
          boot_i <- as.data.frame(rep(i, 4))
          # rename boot column
          names(boot_i) <- "boot_i"
          # make row of all of this info for each bootstrap step
          datBoot_i <- cbind(Time_i, Filter_i, Treatment_i, Size_i, FmI, boot_i, Rep_i)
          
          # Put data in dataframe for later use
          datBoot <- rbind(datBoot, datBoot_i)
        }
      }
    }
  }
}

# create df of zoop density treatments for joining to bootstrapped data
zpdensbms <- sampinfo %>% 
  left_join(ind_dw, by = "Time") %>% 
  # mutate(bms_ugL = Density_numL * ind_dw_ug) %>% 
  select(-c(InitialTime, FinalTime, Volume_mL)) %>% 
  rename(Time_i = Time,
         Treatment_i = Treatment,
         Rep_i = Trial,
         Size_i = Size) %>% 
  mutate(ZP_num_mL = Density_numL / 1000, # convert zp density to mL units (puts CR in terms of mL / ind / d)
         bms_ugCL = Density_numL * Ind_wt_ugC, # get zoop biomass in each treatment
         bms_ugCmL = bms_ugCL / 1000) # convert from L to mL (puts CR in terms of mL / ug zoop C / d)
  
# clean up datBoot df, join zpdens df to calculate CR
datBoot2 <- as.data.frame(datBoot[-1,]) %>% 
  mutate(Time_i = as.factor(Time_i),
         Filter_i= as.factor(Filter_i),
         Treatment_i = as.numeric(Treatment_i),
         Size_i = as.factor(Size_i)) %>% 
  left_join(zpdensbms, by = c("Time_i", "Size_i", "Treatment_i", "Rep_i"))

# fit linear models to estimate slope (clearance rate)
datBootlms <- datBoot2 %>% 
  # create ID variable
  mutate(ID = paste0(Time_i, "_", Size_i, "_", Filter_i, "_", boot_i)) %>% 
  # create list nested by ID
  nest(data = c(-ID)) %>% 
  mutate(
    # fit linear model for each value of ID
    fit = map(data, ~lm(FmI ~ bms_ugCmL, data = .x)), # using ugC/mL for zp bms here b/c we want the CR in terms of mL filtered / ug C zoop / day
    # generate a summary of the linear model summary statistics
    tidied = map(fit, tidy)) 

# extract slopes and intercepts
datBootlms_coefs <- datBootlms %>% 
  # pull info out of tidied component of this list
  unnest(tidied) %>% 
  # rename stuff in term column
  mutate(term = ifelse(term == "(Intercept)","Int", "Slp")) %>% 
  # transform int and slp to wide format
  pivot_wider(id_cols = c(ID, data, fit), names_from = term, values_from = c(estimate:p.value)) %>% 
  # separate ID
  separate(ID, into = c("Time", "Size", "Filter", "BootNum"), sep = "_") %>% 
  # clean up
  select(Time, Size, Filter, BootNum, estimate_Int:p.value_Slp)

# summarize slope and intercept estimates to obtain indiv clearance rates (slope- mL / ugC / day) and max phyto growth rates (intercept- 1 / day) w/ 80% confidence interval
CR_growrate_est <- datBootlms_coefs %>% 
  group_by(Time, Size, Filter) %>% 
  summarize(Int_med = median(estimate_Int, na.rm = TRUE),
            Int_10 = quantile(estimate_Int, 0.1, na.rm = TRUE),
            Int_90 = quantile(estimate_Int, 0.9, na.rm = TRUE),
            CR_mL_ugC_d_med = median(-estimate_Slp, na.rm = TRUE),
            CR_mL_ugC_d_10 = quantile(-estimate_Slp, 0.1, na.rm = TRUE),
            CR_mL_ugC_d_90 = quantile(-estimate_Slp, 0.9, na.rm = TRUE),
            CR_L_ugC_d_med = CR_mL_ugC_d_med / 1000, # convert CR to L units
            CR_L_ugC_d_10 = CR_mL_ugC_d_10 / 1000,
            CR_L_ugC_d_90 = CR_mL_ugC_d_90 / 1000)

# output data
write.csv(CR_growrate_est, "02_CR_growrate_bms_2021_36873_ZPgrazing_boot100.csv")
