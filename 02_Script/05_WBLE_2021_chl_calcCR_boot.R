# Collis/Hood, June 2023
# Updated June 2024 to change units to zooplankton biomass (carbon). Also changing C:chl ratio to 30
# Calculate bootstrapped estimates of indiv zoop clearance rates for chl-a (extracted)
# Doing only for mesozoops- Early Aug and Early Sept for the nano/pico grazing MS

# Libraries ----
library(tidyverse)
library(broom)


# load data ----
# chl data calculated in R script "01_2021_WBLE_ZP_grazingZNR_chla_datamunge"
ZP_chl <- read.csv("01_2021_WBLE_ZP_grazingZNR_chl_wide.csv", row.names = 1) %>% 
  mutate(SampleDate = as.character(as.POSIXct(SampleDate, format = "%Y-%m-%d")))

# selecting only mesozoop Early Aug and Early Sept data (i.e., experiments for nano/pico MS)
# cleaning up df for loop
ZP_chl2 <- ZP_chl %>% 
  filter(Type == "MESO") %>% 
  filter(Event == "Early August 2021" | Event == "Early September 2021") %>% 
  mutate(Time = case_when(Event == "Early August 2021" ~ "EarlyAug",
                          Event == "Early September 2021" ~ "EarlySep",
                          TRUE ~ "blah"),
         Trial = case_when(Rep == 1 ~ "R1",
                           Rep == 2 ~ "R2",
                           Rep == 3 ~ "R3",
                           Rep == 4 ~ "R4",
                           TRUE ~ "blah")) %>% 
  select(Time, Treatment, Trial, Initial_chl_ugL, Final_chl_ugL) %>% 
  rename(Initial = Initial_chl_ugL,
         Final = Final_chl_ugL)

# load sample info from experiments
sampinfo <- read.csv("00_ChlFluoro_WBLE_2021_SampleInfo.csv")

# load field dens and bms-- will use to calculate avg ind dw and apply to zooplankton experiment treatments
# zp_field_bms <- read.csv("Output/csv/00_LEPAS_WBLE_2021_36873_TotBMS_ugCL.csv") %>% 
zp_field_bms <- read.csv("00_LEPAS_WBLE_2021_36873_FieldExp_BMS_ugCL.csv") %>% 
  select(-X)

# calculate avg individual dry weight of zooplankton in each experiment, subset for August and September experiments only
ind_dw <- zp_field_bms %>% 
  filter(Time == "EarlyAug" | Time == "EarlySep") %>% 
  select(Time, Ind_wt_ugC)

# join data for bootstrapping
chlj <- ZP_chl2 %>% 
  left_join(sampinfo, by = c("Time", "Treatment", "Trial"))


#### Estimate clearance rates (CR)-- bootstrap CR estimates by randomly selecting 1 of the 4 Initial/Final bottle reps from each treatment ####
## Loop: calculate r (ln(Final/Initial)/Time) for a given Event/ZP Treatment-- bottle reps randomly generated
dat <- chlj

## create empty dataframe for bootstrapped data ----
datBoot <- as.data.frame(matrix(nrow = 1, ncol = 5))
names(datBoot) <- c("Time_i", "Treatment_i", "FmI", "boot_i", "Rep_i")

# Time_i, Filter_i, Treatment_i, Size_i, FmI, boot_i, Rep_i

## Bootstrapping
# Bootstrap 200 iterations
BootReps <- 100

# select event
for(e in 1:length(unique(dat$Time))){
  # Time_e <- unique(dat$Time)[1]
  Time_e <- unique(dat$Time)[e]
  dat_e <- dat %>% 
    filter(Time == Time_e)
   
      # select treatment
      for(t in 1:length(unique(dat_e$Treatment))){
        # Treatment_t <- unique(dat_e$Treatment)[1]
        Treatment_t <- unique(dat_e$Treatment)[t]
        dat_t <- dat_e %>% 
          filter(Treatment == Treatment_t) %>% 
          mutate(FinalTime = as.POSIXct(FinalTime, format = "%m/%d/%Y %H:%M"),
                 InitialTime = as.POSIXct(InitialTime, format = "%m/%d/%Y %H:%M"))
        
        # within each filter-treatment combine all four finals with randomly selected initials
        for(i in 1:BootReps){
          # i = 1
          # log (Final/initial) where initials are randomly selected from four options within filter-treatment with replacement
          NumSeq_i <- sample(seq(1, 4, by = 1), size = 4, replace = 4)
          
          FmI <- log((dat_t$Final/dat_t[c(NumSeq_i),]$Initial)) / ((as.numeric(dat_t$FinalTime) - as.numeric(dat_t[c(NumSeq_i),]$InitialTime))) * 60 * 60 * 24
          
          # vector of event
          Time_i <- as.data.frame(rep(Time_e, 4))
          # rename event column
          names(Time_i) <- "Time_i"
          # vector of rep numbers
          Rep_i <- dat_t$Trial
          # vector of treatment values
          Treatment_i <- as.data.frame(rep(Treatment_t, 4))
          # rename treatment column
          names(Treatment_i) <- "Treatment_i"
          # boot strap number, ranges from 1 to BootReps
          boot_i <- as.data.frame(rep(i, 4))
          # rename boot column
          names(boot_i) <- "boot_i"
          # make row of all of this info for each bootstrap step
          datBoot_i <- cbind(Time_i, Treatment_i, FmI, boot_i, Rep_i)
          
          # Put data in dataframe for later use
          datBoot <- rbind(datBoot, datBoot_i)
        }
      }
    }

# create df of zoop treatments (by biomass-- carbon) for joining to bootstrapped data
zpbms_exp <- sampinfo %>% 
  select(-c(InitialTime, FinalTime)) %>% 
  left_join(ind_dw, by = "Time") %>% 
  rename(Time_i = Time,
         Treatment_i = Treatment,
         Rep_i = Trial) %>% 
  mutate(ZP_num_mL = Density_numL / 1000, # convert zp density to mL units (puts CR in terms of mL / ind / d)
         bms_ugCL = Density_numL * Ind_wt_ugC, # get zoop biomass in each treatment
         bms_ugCmL = bms_ugCL / 1000) # convert from L to mL (puts CR in terms of mL / ug zoop C / d)

# clean up datBoot df, join zpdens df to calculate CR
datBoot2 <- as.data.frame(datBoot[-1,]) %>% 
  mutate(Time_i = as.factor(Time_i),
         Treatment_i = as.numeric(Treatment_i)) %>% 
  left_join(zpbms_exp, by = c("Time_i", "Treatment_i", "Rep_i"))

# fit linear models to estimate slope (clearance rate)
datBootlms <- datBoot2 %>% 
  # create ID variable
  mutate(ID = paste0(Time_i, "_", boot_i)) %>% 
  # create list nested by ID
  nest(data = c(-ID)) %>% 
  mutate(
    # fit linear model for each value of ID
    fit = map(data, ~lm(FmI ~ bms_ugCmL, data = .x)), # using ugC/mL for zoop biomass here b/c we want the CR in terms of mL / ugC / day
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
  separate(ID, into = c("Time", "BootNum"), sep = "_") %>% 
  # clean up
  select(Time, BootNum, estimate_Int:p.value_Slp)

# summarize slope and intercept estimates to obtain indiv clearance rates (slope- mL / indiv / day) and max phyto growth rates (intercept- 1 / day) w/ 80% confidence interval
CR_growrate_est <- datBootlms_coefs %>% 
  group_by(Time) %>% 
  summarize(Int_med = median(estimate_Int, na.rm = TRUE),
            Int_10 = quantile(estimate_Int, 0.1, na.rm = TRUE),
            Int_90 = quantile(estimate_Int, 0.9, na.rm = TRUE),
            CR_mL_ugC_d_med = median(-estimate_Slp, na.rm = TRUE),
            CR_mL_ugC_d_10 = quantile(-estimate_Slp, 0.1, na.rm = TRUE),
            CR_mL_ugC_d_90 = quantile(-estimate_Slp, 0.9, na.rm = TRUE),
            CR_L_ugC_d_med = CR_mL_ugC_d_med / 1000, # convert CR to L units
            CR_L_ugC_d_10 = CR_mL_ugC_d_10 / 1000,
            CR_L_ugC_d_90 = CR_mL_ugC_d_90 / 1000)

# output data -- 200 boot reps
write.csv(CR_growrate_est, "05_CR_growrate_est_2021_36873_ZPgrazing_chl_bms_boot100.csv")
