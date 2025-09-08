# Collis, June 2023
# updated June 2024, changed to biomass units (ug C)
# Create CR, GR, %SS, upot grazed, selectivity figures

# Libraries ----
library(tidyverse)
library(broom)
library(ggrepel)
library(patchwork)
library(ggpubr)

# load data ----
dat_fin <- read.csv("09_2021_36873_NanoPicoMS_CRGR_final.csv") %>% 
  select(-X)

# reorder factors
dat_fin2 <- dat_fin %>%
  mutate(FoodType = fct_relevel(FoodType, 
                                "Extracted Chl-a", "Cryptophyta and Bacillariophyta",
                                "Chlorophyta", "Cyanobacteria",
                                "Total Nanoplankton","Autotrophic Nanoplankton", 
                                "Heterotrophic Nanoplankton","Total Picoplankton", 
                                "Autotrophic Picoplankton", "Heterotrophic Picoplankton"),
         Date = fct_relevel(Date, "August", "September"),
         FoodType2 = case_when(FoodType == "Extracted Chl-a" ~ "Chl-a",
                               FoodType == "Cryptophyta and Bacillariophyta" ~ "Crypt + Bac",
                               FoodType == "Chlorophyta" ~ "Chloro",
                               FoodType == "Cyanobacteria" ~ "Cyano",
                               FoodType == "Total Nanoplankton" ~ "Tot Nano",
                               FoodType == "Autotrophic Nanoplankton" ~ "Aut Nano",
                               FoodType == "Heterotrophic Nanoplankton" ~ "Het Nano",
                               FoodType == "Total Picoplankton" ~ "Tot Pico",
                               FoodType == "Autotrophic Picoplankton" ~ "Aut Pico",
                               FoodType == "Heterotrophic Picoplankton" ~ "Het Pico",
                               TRUE ~ "NA"),
         FoodType3 = case_when(FoodType == "Extracted Chl-a" ~ "Chl-a",
                               FoodType == "Cryptophyta and Bacillariophyta" ~ "Cr+Ba",
                               FoodType == "Chlorophyta" ~ "Chlor",
                               FoodType == "Cyanobacteria" ~ "Cyano",
                               FoodType == "Total Nanoplankton" ~ "TNano",
                               FoodType == "Autotrophic Nanoplankton" ~ "ANano",
                               FoodType == "Heterotrophic Nanoplankton" ~ "HNano",
                               FoodType == "Total Picoplankton" ~ "TPico",
                               FoodType == "Autotrophic Picoplankton" ~ "APico",
                               FoodType == "Heterotrophic Picoplankton" ~ "HPico",
                               TRUE ~ "NA"),
         Facet2 = case_when(Facet == "Chlorophyll-a" ~ "Total Phytoplankton",
                            Facet == "Fluoroprobe" ~ "FluoroProbe",
                            TRUE ~ Facet),
         Facet2 = fct_relevel(Facet2, "Total Phytoplankton", "FluoroProbe", "Nanoplankton", "Picoplankton"),
         FoodType4 = case_when(FoodType == "Cryptophyta and Bacillariophyta" ~ "Cryptophyta/Bacillariophyta",
                               TRUE ~ FoodType),
         FoodType4 = fct_relevel(FoodType4, 
                                "Extracted Chl-a", "Cyanobacteria",
                                "Cryptophyta/Bacillariophyta",
                                "Chlorophyta", 
                                "Total Nanoplankton","Autotrophic Nanoplankton", 
                                "Heterotrophic Nanoplankton","Total Picoplankton", 
                                "Autotrophic Picoplankton", "Heterotrophic Picoplankton"),
         Date = fct_relevel(Date, "August", "September")) %>% 
  mutate(FoodIDpos = case_when(CR_Ratio_Med > 0 ~ FoodType3,
                            TRUE ~ ""),
         FoodIDneg = case_when(CR_Ratio_Med <= 0 ~ FoodType3,
                               TRUE ~ ""))

#### Make Figures ####

# Figure 4a
# Individual Clearance Rate (mL / ug C/ day)
p1 <-  
ggplot(dat_fin2, aes(x = FoodType4, y = CR_mL_ugC_d_med, group = Date, color = Date)) +        
  geom_point(size = 3, position = position_dodge(width = 0.55),) +
  geom_errorbar(aes(ymin = CR_mL_ugC_d_10, ymax = CR_mL_ugC_d_90), position = position_dodge(width = 0.55), linewidth = 1.1, width = 0.5) +
  facet_grid(.~Facet2, scales = 'free_x')+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  ylab(expression(atop("Clearance rate",
                       paste("(mL µg  ",C^-1," ", d^-1,")")))) +
  ylim(-.3, 12) +
  theme_bw() +
  theme(text = element_text(family = "serif"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 17, color = "black"),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(), 
    legend.text = element_text(size = 18),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 18),
    panel.spacing = unit(0,'lines'))

# Figure 4b
# Community Grazing Rate
p2 <-  
  ggplot(dat_fin2, aes(x = FoodType4, y = GR_ugC_L_d_med, group = Date, color = Date)) +        
  geom_point(size = 3, position = position_dodge(width = 0.55),) +
  geom_errorbar(aes(ymin = GR_ugC_L_d_10, ymax = GR_ugC_L_d_90), position = position_dodge(width = 0.55), linewidth = 1.1, width = 0.5) +
  facet_grid(.~Facet2, scales = 'free_x')+
  ylim(-2, 200)+
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  ylab(expression(atop("Grazing rate",
                       paste("(µg C  ",L^-1," ", d^-1,")"))))  +
  theme_bw() +
  geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  theme(text = element_text(family = "serif"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0,'lines'),
        strip.text = element_blank())

# Figure 4c
# Percent Ambient Biomass Grazed
p3 <- 
  ggplot(dat_fin2, aes(x = FoodType4, y = per_SS_grazed_med, group = Date, color = Date)) +        
  geom_point(size = 3, position = position_dodge(width = 0.55),) +
  geom_errorbar(aes(ymin = per_SS_grazed_10, ymax = per_SS_grazed_90), position = position_dodge(width = 0.55), linewidth = 1.1, width = 0.5) +
  ylim(-3, 80)+
  geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  facet_grid(.~Facet2, scales = 'free_x') +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  ylab(expression(atop("% biomass grazed"))) +
  theme_bw() +
  theme(text = element_text(family = "serif"),
        axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        legend.title= element_blank(), 
        legend.text=element_text(size = 17),
        strip.text = element_blank(),
        panel.spacing = unit(0,'lines'))

#Combining individual grazing and population grazing into one image and then saving
p123 <- p1 / p2 / p3
# 
ggsave('10_WBLE_2021_Fig4_bms_083124.png', p123,
        width = 10, height = 10)


# Figure 5a
# Selectivity graph-- Relative Food Abundance vs. Relative Clearance Rate
p5a <- ggplot(dat_fin2, aes(x = Foodi_FoodMax, y = CR_Ratio_Med, color = Date, label = FoodType3, family = "serif", axes = FALSE,
                            main = "axes = FALSE")) + #size=Biomass
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1, colour="#000100"), linewidth = 1)+
  labs(title=" ")+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0, colour="#000100"), linewidth = 1)+
  labs(title=" ")+
  geom_segment(aes(x = 1, y = 1, xend = 0, yend = 1, colour="#000100"), linewidth = 1)+
  labs(title=" ")+
  geom_segment(aes(x = 1, y = 1, xend = 1, yend = 0, colour="#000100"), linewidth = 1)+
  labs(title=" ")+
  geom_point(size=3) +
  scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  scale_color_manual(values = c("grey20", "#E69F00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  geom_label_repel(size = 6, max.overlaps=Inf,min.segment.length = 0, fill = "white", fontface = "bold", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))+
  xlim(-1, 1.2)+
  ylim(-1, 1.2)+
  ylab("Relative clearance rate")+
  xlab("Relative food abundance")+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1, colour="#000100"),linewidth=1, linetype="dashed")+
  theme_bw()+
  theme(text = element_text(family = "serif"),
    axis.title.y = element_text(size = 14),
    axis.line.y = element_blank(),
    axis.line.x.top = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none",
    legend.text = element_text(size = 12))

ggsave('10_WBLE_2021_Fig5a_062424.png', p5a,
       width = 16, height = 16)

# Figure 5b
# Food "Importance" vs. "Selectivity"-- Relative Clearance Rate vs. Relative Grazing rate
p5b <- ggplot(dat_fin2, aes(x = GR_Ratio_Med, y = CR_Ratio_Med, color = Date, label = FoodType3, family = "serif", size = 5, axes = FALSE,
                            main = "axes = FALSE")) + #size=Biomass
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1, colour="#000100"), linewidth = 1)+
  labs(title=" ")+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0, colour="#000100"), linewidth = 1)+
  labs(title=" ")+
  geom_segment(aes(x = 1, y = 1, xend = 0, yend = 1, colour="#000100"),linewidth = 1)+
  labs(title=" ")+
  geom_segment(aes(x = 1, y = 1, xend = 1, yend = 0, colour="#000100"),linewidth = 1)+
  labs(title=" ")+
  geom_point(size = 3) +
  coord_cartesian(clip = "off") +
  geom_label_repel(size = 6, max.overlaps=Inf,min.segment.length = 0, fill = "white", fontface = "bold", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))+
  scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  scale_color_manual(values = c("grey20", "#E69F00", "#0072B2"), guide = "none") +
  xlim(-1, 1.2)+
  ylim(-1, 1.2)+
  ylab("Relative clearance rate")+
  xlab("Relative grazing rate")+
  geom_segment(aes(x = .5, y = 0, xend = .5, yend = 1, colour="#000100"),linewidth = 1, linetype="dashed")+
  geom_segment(aes(x = 0, y = .5, xend = 1, yend = .5, colour="#000100"),linewidth = 1, linetype="dashed")+
  theme_bw()+
  theme(
    axis.title.y = element_text(size = 14),
    axis.line.y = element_blank(),
    axis.line.x.top = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none",
    legend.text = element_text(size = 12))

ggsave('10_WBLE_2021_Fig5b_062424.png', p5b,
       width = 16, height = 16)

## NEW Figure 6
# Population CR (grazing mortality) vs. Food Growth Rate Figure
# Population CR is the ind CR x Zoop biomass at the field site
# Food Growth Rate is the y-Intercept of the Zoop Treatment-Food Change curve
# Originally, we calculated grazing mortality using the avg ZP biomass in the experiment, but now we're using ZP field site bms
p6 <- ggplot(dat_fin2, aes(x = Int_med, y = m_fld_med, shape = Date, color = Date, family = "serif", label = FoodType3)) +
  geom_errorbar(aes(ymin = m_fld_10, ymax = m_fld_90)) +
  geom_errorbarh(aes(xmin = Int_10, xmax = Int_90)) +
  geom_hline(yintercept = 0, linewidth = 1) +
  geom_vline(xintercept = 0, linewidth = 1) +
  geom_segment(x = 0, y = 0, xend = 1.1, yend = 1.1,
               col = "dark grey", linewidth = 1, linetype = "dashed") + # 100% standing stock grazed
  geom_segment(x = 0, y = 0, xend = 1.0, yend = 1.5,
               col = "dark grey", linewidth = 1, linetype = "dashed") +
  geom_segment(x = 0, y = 0, xend = 2.0, yend = 1.0,
               col = "dark grey", linewidth = 1, linetype = "dashed") +
  geom_label(aes(size = 15, fill = factor(Date)), color = "white", fontface = "bold", position=position_jitter()) +
  scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  labs( x = expression(paste("Food growth rate  ", (d^{-1}))), 
        y = expression(paste("Grazing mortality rate  ", (d^{-1})))) +
  xlim(-0.66, 1.04) +
  ylim(-0.7, 0.75) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  annotate(geom="text", x = 0.65, y = 0.57, label = "100%",
           color = "dark grey", size = 7, family = "serif") +
  annotate(geom="text", x = 0.85, y = 0.38, label = "50%",
           color = "dark grey", size = 7, family = "serif") +
  annotate(geom="text", x = 0.5, y = 0.64, label = "150%",
           color = "dark grey", size = 7, family = "serif") +
  annotate(geom="text", x = .95, y = 0.68, label = "% growth \n consumed",
           color = "black", size = 7,
           fontface = "bold.italic", family = "serif") +
  annotate(geom="text", x = -0.2, y = 0.66, label = "Positive grazing, \n food decline",
           color = "black", size = 7,
           fontface = "bold.italic", family = "serif") +
  annotate(geom="text", x = -0.19, y = -0.6, label = "No grazing, \n food decline",
           color = "black", size = 7,
           fontface = "bold.italic", family = "serif") +
  annotate(geom="text", x = 0.18, y = -0.6, label = "No grazing, \n food increase",
           color = "black", size = 7,
           fontface = "bold.italic", family = "serif") +
  theme_bw() +
  theme(text = element_text(family = "serif"),
        axis.text.x = element_text(size = 23),
        axis.title.x = element_text(size = 23),
        axis.text.y = element_text(size=23),
        axis.title.y = element_text(size=23),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank(), 
        legend.text = element_text(size = 22),
        panel.spacing = unit(0,'lines'),
        strip.text = element_blank(),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.13, .95))

p6
ggsave('10_WBLE_2021_Fig6_new.png', p6,
       width = 14, height = 8)


## OLD Figure 6
# Population CR vs. Food Growth Rate Figure
# Population CR is the ind CR x Avg Zoop biomass in Experiment
# Food Growth Rate is the y-Intercept of the Zoop Treatment-Food Change curve
p6 <- ggplot(dat_fin2, aes(x = Int_med, y = m_med, shape = Date, color = Date, family = "serif", label = FoodType3)) +
      geom_errorbar(aes(ymin = m_10, ymax = m_90)) +
      geom_errorbarh(aes(xmin = Int_10, xmax = Int_90)) +
      geom_hline(yintercept = 0, linewidth = 1) +
      geom_vline(xintercept = 0, linewidth = 1) +
      geom_segment(x = 0, y = 0, xend = 1.1, yend = 1.1,
                   col = "dark grey", linewidth = 1, linetype = "dashed") + # 100% standing stock grazed
      geom_segment(x = 0, y = 0, xend = 1.0, yend = 1.5,
                   col = "dark grey", linewidth = 1, linetype = "dashed") +
      geom_segment(x = 0, y = 0, xend = 2.0, yend = 1.0,
                   col = "dark grey", linewidth = 1, linetype = "dashed") +
      geom_label(aes(size = 15, fill = factor(Date)), color = "white", fontface = "bold", position=position_jitter()) +
      scale_fill_manual(values = c("#E69F00", "#0072B2"), guide = "none") +
      scale_color_manual(values = c("#E69F00", "#0072B2")) +
      labs( x = expression(paste("Food growth rate  ", (d^{-1}))), 
            y = expression(paste("Grazing mortality rate  ", (d^{-1})))) +
      xlim(-0.67, 1.04) +
      ylim(-0.6, 0.67) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      annotate(geom="text", x = 0.65, y = 0.57, label = "100%",
               color = "dark grey", size = 7, family = "serif") +
      annotate(geom="text", x = 0.85, y = 0.38, label = "50%",
               color = "dark grey", size = 7, family = "serif") +
      annotate(geom="text", x = 0.5, y = 0.64, label = "150%",
               color = "dark grey", size = 7, family = "serif") +
      annotate(geom="text", x = .95, y = 0.62, label = "% growth \n consumed",
           color = "black", size = 7,
           fontface = "bold.italic", family = "serif") +
     annotate(geom="text", x = -0.2, y = 0.66, label = "Positive grazing, \n food decline",
           color = "black", size = 7,
           fontface = "bold.italic", family = "serif") +
     annotate(geom="text", x = -0.19, y = -0.57, label = "No grazing, \n food decline",
           color = "black", size = 7,
           fontface = "bold.italic", family = "serif") +
    annotate(geom="text", x = 0.18, y = -0.57, label = "No grazing, \n food increase",
             color = "black", size = 7,
             fontface = "bold.italic", family = "serif") +
        theme_bw() +
      theme(text = element_text(family = "serif"),
            axis.text.x = element_text(size = 23),
            axis.title.x = element_text(size = 23),
            axis.text.y = element_text(size=23),
            axis.title.y = element_text(size=23),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "top",
            legend.title = element_blank(), 
            legend.text = element_text(size = 22),
            panel.spacing = unit(0,'lines'),
            strip.text = element_blank(),
            plot.tag = element_text(size = 16),
            plot.tag.position = c(0.13, .95))

p6
ggsave('10_WBLE_2021_Fig6.png', p6,
       width = 14, height = 8)
