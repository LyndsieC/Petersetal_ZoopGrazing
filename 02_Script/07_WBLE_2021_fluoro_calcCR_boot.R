# Collis, June 2023
# Calculate bootstrapped estimates of indiv zoop clearance rates for fluoroprobe phytoplankton groups
# Doing only for mesozoops- Early Aug and Early Sept for the nano/pico grazing MS

# Libraries ----
library(tidyverse)
library(broom)


# load data ----
# fluoroprobe data
ZP_fluoro <- read.csv("03_fluoroprobe_exp_phytogrps_WBLE2021.csv", row.names = 1) %>% 
  mutate(SampleDate = as.character(as.POSIXct(SampleDate, format = "%Y-%m-%d")))

# selecting only mesozoop Early Aug and Early Sept data (i.e., experiments for nano/pico MS)
# cleaning up df for loop
ZP_fluoro2 <- ZP_fluoro %>% 
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
  select(Time, Treatment, Trial, Phyto_grp, Initial_chl_ugL, Final_chl_ugL) %>% 
  rename(Initial = Initial_chl_ugL,
         Final = Final_chl_ugL)

# Make dataframe wide so I can combine diatoms and cryptophytes
ZP_fluorow <- ZP_fluoro2 %>% 
  pivot_wider(names_from = Phyto_grp, values_from = c(Initial, Final))

# Combine diatom and cryptophytes into one group
ZP_fluorow2 <- ZP_fluorow %>% 
  mutate(Initial_CryptoDiatoms = Initial_Cryptophyta + Initial_Diatoms,
         Final_CryptoDiatoms = Final_Cryptophyta + Final_Diatoms) %>% 
  select(-c(Initial_Cryptophyta, Initial_Diatoms, Final_Cryptophyta, Final_Diatoms)) # remove the individual crypto/diatom groups

# Make data long again
ZP_fluorol <- ZP_fluorow2 %>% 
  pivot_longer(cols = Initial_Bluegreen:Final_CryptoDiatoms,
               names_to = c("In_Fin", "Phyto_grp"),
               names_sep = "_",
               values_to = "chl_ugL")

# pivot wider so I can run bootstrap function
ZP_fluorow3 <- ZP_fluorol %>% 
  pivot_wider(names_from = In_Fin, values_from = chl_ugL)

# load sample info
sampinfo <- read.csv("00_ChlFluoro_WBLE_2021_SampleInfo.csv")

# load field dens and bms-- will use to calculate avg ind dw and apply to zooplankton experiment treatments
zp_field_bms <- read.csv("00_LEPAS_WBLE_2021_36873_FieldExp_BMS_ugCL.csv") %>% 
  select(-X)

# calculate avg individual dry weight of zooplankton in each experiment, subset for August and September experiments only
ind_dw <- zp_field_bms %>% 
  filter(Time == "EarlyAug" | Time == "EarlySep") %>% 
  select(Time, Ind_wt_ugC)

# join fluoroprobe and sampinfo data for bootstrapping
ZP_fluoroj <- ZP_fluorow3 %>% 
  left_join(sampinfo, by = c("Time", "Treatment", "Trial"))

#### Estimate clearance rates (CR)-- bootstrap CR estimates by randomly selecting 1 of the 4 Initial/Final bottle reps from each treatment ####
## Loop: calculate r (ln(Final/Initial)/Time) for a given Event/ZP Treatment/Phyto group-- bottle reps randomly generated
dat <- ZP_fluoroj

## create empty dataframe for bootstrapped data ----
FluoroBoot <- as.data.frame(matrix(nrow = 1, ncol = 6))
names(FluoroBoot) <- c("Time_i", "Treatment_i", "Phyto_grp_i", "FmI", "boot_i", "Rep_i")

# Bootstrapping reps
BootReps <- 100

# bootstraping for loops

# select sample date
for(f in 1:length(unique(dat$Time))){
  # f = 1; z = 1; t = 2; p = 1
  # SampleDate_f <-  unique(Fluoroj$SampleDate)[1]
  Time_f <-  unique(dat$Time)[f]
  dat_f <- dat %>% 
    filter(Time == Time_f)
  
  # select PP group
  for(p in 1:length(unique(dat_f$Phyto_grp))) {
    Phyto_grp_fp <- unique(dat_f$Phyto_grp)[p]
    dat_fp <- dat_f %>% 
      filter(Phyto_grp == Phyto_grp_fp)
    
    # select Treatment
    for(t in 1:length(unique(dat_fp$Treatment))){
      # Treatment_t <- unique(Fluoro$Treatment)[1]
      Treatment_t <- unique(dat_fp$Treatment)[t]
      dat_fpt <- dat_fp %>% 
        filter(Treatment == Treatment_t) %>% 
        mutate(FinalTime = as.POSIXct(FinalTime, format = "%m/%d/%Y %H:%M"),
               InitialTime = as.POSIXct(InitialTime, format = "%m/%d/%Y %H:%M"))
      
      # within each filter-treatment combine all four finals with randomly selected initials
      for(i in 1:BootReps){
        # i = 1
        # log (Final/initial) where initials are randomly selected from four options within filter-treatment with replacement
        NumSeq_i <- sample(seq(1, 4, by = 1), size = 4, replace = 4)
        
        FmI <- log((dat_fpt$Final/dat_fpt[c(NumSeq_i),]$Initial)) / ((as.numeric(dat_fpt$FinalTime) - as.numeric(dat_fpt[c(NumSeq_i),]$InitialTime))) * 60 * 60 * 24
        
        # vector of event
        Time_i <- as.data.frame(rep(Time_f, 4))
        # rename event column
        names(Time_i) <- "Time_i"
        # vector of rep numbers
        Rep_i <- dat_fpt$Trial
        # vector of treatment values
        Treatment_i <- as.data.frame(rep(Treatment_t, 4))
        # rename treatment column
        names(Treatment_i) <- "Treatment_i"
        # vector of phyto groups
        Phyto_grp_i <- as.data.frame(rep(Phyto_grp_fp,4))
        # rename phyto group
        names(Phyto_grp_i) <- "Phyto_grp_i"
        # boot strap number, ranges from 1 to BootReps
        boot_i <- as.data.frame(rep(i, 4))
        # rename boot column
        names(boot_i) <- "boot_i"
        # make row of all of this info for each bootstrap step
        FluoroBoot_i <- cbind(Time_i, Treatment_i, Phyto_grp_i, FmI, boot_i, Rep_i)
        
        # Put data in dataframe for later use
        FluoroBoot <- rbind(FluoroBoot, FluoroBoot_i)
      }
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
FluoroBoot2 <- as.data.frame(FluoroBoot[-1,]) %>% 
  mutate(Time_i = as.factor(Time_i),
         Phyto_grp_i = as.factor(Phyto_grp_i),
         Treatment_i = as.numeric(Treatment_i)) %>% 
  left_join(zpbms_exp, by = c("Time_i", "Treatment_i", "Rep_i"))

# fit linear models to estimate slope (clearance rate)
FluoroBootlms <- FluoroBoot2 %>% 
  # create ID variable
  mutate(ID = paste0(Time_i, "_", Phyto_grp_i, "_", boot_i)) %>% 
  # create list nested by ID
  nest(data = c(-ID)) %>% 
  mutate(
    # fit linear model for each value of ID
    fit = map(data, ~lm(FmI ~ bms_ugCmL, data = .x)), # using ugC/mL for zoop biomass here b/c we want the CR in terms of mL / ugC / day
    # generate a summary of the linear model summary statistics
    tidied = map(fit, tidy)) 

# extract slopes and intercepts
FluoroBootlms_coefs <- FluoroBootlms %>% 
  # pull info out of tidied component of this list
  unnest(tidied) %>% 
  # rename stuff in term column
  mutate(term = ifelse(term == "(Intercept)","Int", "Slp")) %>% 
  # transform int and slp to wide format
  pivot_wider(id_cols = c(ID, data, fit), names_from = term, values_from = c(estimate:p.value)) %>% 
  # separate ID
  separate(ID, into = c("Time", "Phyto_Grp", "BootNum"), sep = "_") %>% 
  # clean up
  select(Time, Phyto_Grp, BootNum, estimate_Int:p.value_Slp)

# summarize slope and intercept estimates to obtain indiv clearance rates (slope- mL / indiv / day) and max phyto growth rates (intercept- 1 / day) w/ 80% confidence interval
fluoro_CR_growrate_est <- FluoroBootlms_coefs %>% 
  group_by(Time, Phyto_Grp) %>% 
  summarize(Int_med = median(estimate_Int, na.rm = TRUE),
            Int_10 = quantile(estimate_Int, 0.1, na.rm = TRUE),
            Int_90 = quantile(estimate_Int, 0.9, na.rm = TRUE),
            CR_mL_ugC_d_med = median(-estimate_Slp, na.rm = TRUE),
            CR_mL_ugC_d_10 = quantile(-estimate_Slp, 0.1, na.rm = TRUE),
            CR_mL_ugC_d_90 = quantile(-estimate_Slp, 0.9, na.rm = TRUE),
            CR_L_ugC_d_med = CR_mL_ugC_d_med / 1000, # convert CR to L units
            CR_L_ugC_d_10 = CR_mL_ugC_d_10 / 1000,
            CR_L_ugC_d_90 = CR_mL_ugC_d_90 / 1000)

# output data-- 100 boot reps
write.csv(fluoro_CR_growrate_est, "07_CR_growrate_est_2021_36873_ZPgrazing_fluoro_boot100.csv")
