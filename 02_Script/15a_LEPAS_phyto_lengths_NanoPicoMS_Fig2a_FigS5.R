# load libraries
library(tidyverse)
library(ggdist)
library(gghalves)
library(RColorBrewer)

# read in data
# density/bms for each phyto taxa
phyto <- read.csv("03b_LEPAS_phyto_biomass_09012024.csv") %>% 
  select(-X) %>% 
  mutate(Sample_date = as.POSIXct(Sample_date, format = "%m/%d/%Y"),
         Type = as.factor(Type),
         Edible = as.factor(Edible),
         Size_class = as.factor(Size_class),
         Length2 = replace_na(Length_um, 0),
         Diameter2 = replace_na(Diameter_cell_um, 0),
         Length_um3 = Length2 + Diameter2)
str(phyto)

# calc total per phylum
phyto_tot <- phyto %>% 
  group_by(Sample_date, Phylum) %>% 
  summarize(bms_tot = sum(Phyto_biomass_mgL, na.rm = TRUE),
            dens_tot = sum(Density_cellsPerL, na.rm = TRUE))

# calculate percent genus in each phylum
per_genus <- phyto %>% 
  left_join(phyto_tot, by = c("Sample_date", "Phylum")) %>% 
  mutate(per_dens = (Density_cellsPerL / dens_tot) * 100,
         per_bms = (Phyto_biomass_mgL / bms_tot) * 100)

# pull out just cryptos and diatoms, get total density and bms
cryp_dia_tot <- phyto %>% 
  filter(Phylum == "Bacillariophyta" | Phylum == "Cryptophyta") %>% 
  group_by(Sample_date) %>% 
  summarize(bms_tot = sum(Phyto_biomass_mgL, na.rm = TRUE),
            dens_tot = sum(Density_cellsPerL, na.rm = TRUE))

# get per Aula of the Cryptos/Diatoms
per_genus_CD <- phyto %>% 
  filter(Phylum == "Bacillariophyta" | Phylum == "Cryptophyta") %>% 
  left_join(cryp_dia_tot, by = c("Sample_date")) %>% 
  mutate(per_dens = (Density_cellsPerL / dens_tot) * 100,
         per_bms = (Phyto_biomass_mgL / bms_tot) * 100)

# individual phyto length data
# we need this because some phyto are >200um-- these need to be removed when comparing to fluoroprobe data from the experiment (water was sieved through a 200um sieve)
phyto_ind <- read.csv("03b_LEPAS_phyto_RawLengths.csv") %>%
  select(-X) %>% 
  mutate(SampleDate = as.POSIXct(SampleDate, format = "%m/%d/%Y"),
         Ct = 1)

#### Individual Phyto Lengths ####
# Need to account for phyto >200um to make the 200um-sieved fluoroprobe data directly comparable
# Only phyto >200 are Aulacoseira and Chlorophyte filament, so those are the taxa that need to be "corrected"
phyto_over200 <- phyto_ind %>% 
  filter(Size_cat == "Lrg_200")

# Calculate proportion of Aulacoseira above and below 200um
Aulo <- phyto_ind %>% 
  filter(Genus_type == "Aulacoseira") %>%
  filter(Length_um > 0) %>% 
  mutate(BV = (3.14 * Width1_um^2 * Length_um) / 4)

Aulo_szgrp <- Aulo %>% 
  group_by(Sample_ID, Size_cat) %>% 
  summarise(BV = sum(BV, na.rm = TRUE),
            Ct = sum(Ct, na.rm = TRUE))

Aulo_tot <- Aulo %>% 
  group_by(Sample_ID) %>% 
  summarise(BV = sum(BV, na.rm = TRUE),
            Ct = sum(Ct, na.rm = TRUE))

Aulo_per <- left_join(Aulo_szgrp, Aulo_tot, by = "Sample_ID") %>% 
  mutate(perBV = BV.x / BV.y,
         perCt = Ct.x / Ct.y)

# Calculate proportion of Chlorophyte filaments above and below 200um
Chloro <- phyto_ind %>% 
  filter(Genus_type == "Chlorophyte_Filament") %>% 
  filter(Length_um > 0) %>% 
  mutate(BV = (3.14 * Width1_um^2 * Length_um) / 4)

Chloro_szgrp <- Chloro %>% 
  group_by(Sample_ID, Size_cat) %>% 
  summarise(BV = sum(BV, na.rm = TRUE),
            Ct = sum(Ct, na.rm = TRUE))

Chloro_tot <- Chloro %>% 
  group_by(Sample_ID) %>% 
  summarise(BV = sum(BV, na.rm = TRUE),
            Ct = sum(Ct, na.rm = TRUE))

Chloro_per <- left_join(Chloro_szgrp, Chloro_tot, by = "Sample_ID") %>% 
  mutate(perBV = BV.x / BV.y,
         perCt = Ct.x / Ct.y)

# create column to proportion Aulacoseira and Chlorophyte filaments so biomass is "corrected" for only organisms <200um
# "per" is taken from "Aulo_per" and "Chloro_per" (perBV column)
phyto2 <- phyto %>% 
  mutate(per = case_when(Sample_ID == "20210809-36-873" & Genus_type == "Aulacoseira" ~ 0.2649811,
                         Sample_ID == "20210910-36-873" & Genus_type == "Aulacoseira" ~ 0.2665956,
                         Sample_ID == "20210809-36-873" & Genus_type == "Chlorophyte Filament" ~ 0.3950617,
                         TRUE ~ 1), # everything else will remain unchanged b/c all other individuals are <200um
         BMS_mgL_under200 = Phyto_biomass_mgL * per) 

# Subset fluoroprobe groups
phyto2a <- phyto2 %>% 
  filter(Phylum == "Chlorophyta" | Phylum == "Cryptophyta" | Phylum == "Bacillariophyta" | Phylum == "Cyanobacteria")

# summarize bms by phylum for all phyto and only phyto <200um
phyto_phyl_sum <- phyto2a %>% 
  group_by(Sample_date, Phylum) %>% 
  summarize(BMS_mgL_all = sum(Phyto_biomass_mgL, na.rm = TRUE), # this is phyto of all sizes
            BMS_mgL_less200 = sum(BMS_mgL_under200, na.rm = TRUE)) # this is only phyto under 200um

phyto_phyl_sum2 <- phyto_phyl_sum %>% 
  mutate(WW_C_conv = case_when(Phylum == "Bacillariophyta" | Phylum == "Cryptophyta" ~ "0.11", # using diatom WW-C conv factor for diatoms for diatoms and cyptos
                               Phylum == "Chlorophyta" ~ "0.16",
                               Phylum == "Cyanobacteria" ~ "0.22",
                               TRUE ~ "NA")) %>% 
  mutate(bms_ugC_less200_tot = (BMS_mgL_less200 * as.numeric(WW_C_conv)) * 1000)
  
write.csv(phyto_phyl_sum2, "15a_LEPAS_phyto_bms_less200um_2021_36873.csv")

# convert bms in phyto2a to C units
phyto_less200_Cunits <- phyto2a %>% 
  mutate(WW_C_conv = case_when(Phylum == "Bacillariophyta" | Phylum == "Cryptophyta" ~ "0.11", # using diatom WW-C conv factor for diatoms for diatoms and cyptos
                               Phylum == "Chlorophyta" ~ "0.16",
                               Phylum == "Cyanobacteria" ~ "0.22",
                               TRUE ~ "NA")) %>% 
  mutate(bms_ugC_less200 = (BMS_mgL_under200 * as.numeric(WW_C_conv)) * 1000)

# calculate percent of each genus towards each phylum total
phyto_less200_per <- phyto_less200_Cunits %>% 
  left_join(phyto_phyl_sum2 %>% select(c(Sample_date, Phylum, bms_ugC_less200_tot)), by = c("Sample_date", "Phylum")) %>% 
  mutate(bms_per = (bms_ugC_less200 / bms_ugC_less200_tot) * 100,
         Date = case_when(Sample_date == "2021-08-09" ~ "August",
                          Sample_date == "2021-09-10" ~ "September",
                          TRUE ~ "NA"),
         Size_class2 = case_when(Size_class == "MICRO" ~ "Micro",
                                 Size_class == "NANO" ~ "Nano",
                                 Size_class == "PICO" ~ "Pico",
                                 TRUE ~ "NA")) %>% 
  rename(Taxa = Genus_type) %>% 
  mutate(Taxa2 = case_when(Phylum == "Chlorophyta" & bms_per < 5 ~ "Other Chlorophyta",
                           Phylum == "Bacillariophyta" & bms_per < 5 ~ "Other Bacillariophyta",
                           Phylum == "Cryptophyta" & bms_per < 5 ~ "Other Cryptophyta",
                           Phylum == "Cyanobacteria" & bms_per < 5 ~ "Other Cyanobacteria",
                           TRUE ~ Taxa)) %>% 
  mutate(Phylum = factor(Phylum, levels = c("Cyanobacteria", "Cryptophyta", "Bacillariophyta", "Chlorophyta")),
         Taxa2 = factor(Taxa2, levels = c("Microcystis", "Other Cyanobacteria", 
                                        "Chroomonas", "Cryptomonas", "Rhodomonas", "Other Cryptophyta",
                                        "Aulacoseira", "Centric Diatom", "Cyclotella", "Other Bacillariophyta",
                                        "Carteria", "Chlamydomonas", "Chlorophyte Filament", "Eudorina", "Solitary Green", "Spiny Green", "Other Chlorophyta")))

write.csv(phyto_less200_per, "15a_LEPAS_phyto_bms_less200um_per_2021_36873.csv")

# sum all genuses where biomass is <5% as "other"
phyto_less200_per_othersum <- phyto_less200_per %>% 
  group_by(Date, Phylum, Taxa2) %>% 
  summarize(bms_ugC_sum = sum(bms_ugC_less200, na.rm = TRUE))

## Plot biomass of each phytoplankton genus-- all taxa <5% biomass as listed as "Other"
png(file = "15a_LEPASphyto_AugSep_bms_less200um_FigS5.png",
    width = 600, height = 750)

ggplot(phyto_less200_per_othersum, aes(y = bms_ugC_sum, x = Phylum, fill = Taxa2)) +
  geom_bar(stat = "identity", colour="black") +
  ylab(expression(paste("Phytoplankton Biomass (µg C  ",  L^-1,")"))) +
  facet_grid(.~Date, scales = "free") +
  scale_fill_manual(values = c("green4", "grey67", # cyano colors
                               "goldenrod3", "blue2", "brown1",  "black", # crypto colors
                               "purple2", "plum1", "cyan2", "ivory2", # diatom colors
                               "chartreuse", "royalblue1", "orange2", "violetred", "cadetblue2", "springgreen4", "khaki4")) + # Chlorophyta colors
  theme_light() +
  theme(text = element_text(family = "serif"),
        axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(size=15, color = "black"),
        axis.title.y = element_text(size=17, color = "black"),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.title = element_blank())

dev.off()