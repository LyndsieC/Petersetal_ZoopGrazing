# get data from 14 and 15b

phyto_comb <- bind_rows(comb_dat, field_chl) %>% 
  mutate(Type = as.factor(Type),
         Type = fct_relevel(Type, "Chlorophyll-a", "FluoroProbe", "Microscope"),
         Phyto_grp2 = as.factor(Phyto_grp2),
         Phyto_grp2 = fct_relevel(Phyto_grp2, "Total phytoplankton", "Bacillariophyta", "Cryptophyta",
                                  "Cryptophyta and Bacillariophyta", "Chlorophyta", "Cyanobacteria"))
         
# Make plot of all three phyto measurements
ggplot(phyto_comb, aes(y = bms_ugCL, x = Event, fill = Phyto_grp2)) +
  geom_bar(stat = "identity") +
  facet_wrap(vars(Type), scales = "free_y") + # this allows the y-axis to be scaled differently for fluoroprobe vs. microscope
  scale_fill_manual(values = c("Total phytoplankton" = "gray23",
                               "Cyanobacteria" = "deepskyblue3",
                               "Chlorophyta" = "chartreuse3",
                               "Cryptophyta" = "darkgoldenrod1",
                               "Bacillariophyta" = "darkgoldenrod3",
                               "Cryptophyta and Bacillariophyta" = "darkred")) +
  scale_x_discrete(labels = c("August", "September")) +
  ylab(expression(paste("Phytoplankton Biomass (µg C  ", L^-1,")"))) +
  theme_bw() +
  theme(text = element_text(family = "serif"),
        axis.title.y = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 23),
        legend.position = "top",
        legend.title= element_blank(), 
        panel.background = element_rect(fill = "transparent", linewidth = 1),
        axis.ticks = element_line(color = "black"),
        panel.grid = element_line(color = "transparent"),
        strip.text = element_text(size = 25),
        strip.background.y = element_rect(fill = "transparent", color = "transparent", linewidth = 1),
        strip.background.x = element_rect(fill = "transparent"),
        plot.tag = element_text(size = 16),) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
