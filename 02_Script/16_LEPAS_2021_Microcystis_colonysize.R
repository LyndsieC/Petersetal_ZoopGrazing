library(tidyverse)

dat <- read.csv("01_Data/00_Microcystis_ColonySize_Measurements.csv") %>% 
        mutate(Biomass = (pi*Dcol)/6) %>% 
        group_by(Sample_ID, NewOld) %>% 
        mutate(PerBiomass = Biomass/sum(Biomass))

bms_sum <- dat %>% 
  group_by(Sample_ID, NewOld) %>% 
  summarize(MicrocystisBiomassTot = sum(Biomass, na.rm = TRUE)) %>% 
  mutate(Bms_tot = MicrocystisBiomassTot * 0.22)

ggplot(dat, aes(y = Dcol, x = PerBiomass)) +
  geom_point()

ggplot(dat %>% 
         filter(NewOld == "New"), 
       aes(x = Dcol)) +
  geom_histogram(fill = "dodgerblue", color = "black") +
  scale_x_log10() +
  facet_wrap(~Sample_ID) +
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 100)


# % of COLONIES by biomass
# <
L200 <- dat %>% 
  mutate(L200 = ifelse(Dcol < 200, "L200", "G200"),
         L60 = ifelse(Dcol < 60, "L60", "G60")) %>% 
  group_by(Sample_ID, NewOld, L200) %>% 
  summarize(Biomass = sum(Biomass)) %>% 
  pivot_wider(id_cols = Sample_ID:NewOld, names_from = "L200", values_from = Biomass) %>% 
  mutate(LessThan = L200/(L200 + G200),
         Frac = "L200") %>% 
  filter(NewOld == "New") %>% 
  ungroup() %>% 
  select(-NewOld, -G200, -L200)

L100 <- dat %>% 
  mutate(L200 = ifelse(Dcol < 200, "L200", "G200"),
         L100 = ifelse(Dcol < 100, "L100", "G100")) %>% 
  group_by(Sample_ID, NewOld, L100) %>% 
  summarize(Biomass = sum(Biomass)) %>% 
  pivot_wider(id_cols = Sample_ID:NewOld, names_from = "L100", values_from = Biomass) %>% 
  mutate(LessThan = L100/(L100 + G100),
         Frac = "L100")%>% 
  filter(NewOld == "New")%>% 
  ungroup() %>% 
  select(-NewOld, -G100, -L100)


L20 <- dat %>% 
  mutate(L200 = ifelse(Dcol < 200, "L200", "G200"),
         L20 = ifelse(Dcol < 20, "L20", "G20")) %>% 
  group_by(Sample_ID, NewOld, L20) %>% 
  summarize(Biomass = sum(Biomass)) %>% 
  pivot_wider(id_cols = Sample_ID:NewOld, names_from = "L20", values_from = Biomass) %>% 
  mutate(LessThan = L20/(L20 + G20),
         Frac = "L20")%>% 
  filter(NewOld == "New")%>% 
  ungroup() %>% 
  select(-NewOld, -G20, -L20)

MicrocystisFracs <- rbind(L200, L100, L20) %>% 
  pivot_wider(id_cols = Sample_ID, names_from = Frac, values_from = LessThan) %>% 
  mutate(FG200 = 1-L200,
         F100to200 = L200-L100,
         F20to100 = L100-L20) %>% 
  select(-L200, -L100) %>% 
  pivot_longer(cols = L20:F20to100) %>% 
  mutate(MicrocystisBiomassTot = ifelse(Sample_ID == "36-873_08092021", 1.8, 1.14), # 1.8 and 1.14 is the total Microcystis biomass from LEPAS on each date
         MicrocystisBiomassFrac = MicrocystisBiomassTot*value,
         name = fct_relevel(as.factor(name), "FG200", "F100to200", "F20to100", "L20"),
         name = fct_recode(name, 
                           "< 20 µm" = "L20",
                           "20 <- 100 µm" = "F20to100",
                           "100 - 200 µm" = "F100to200",
                           ">200 µm" = "FG200"),
         Month = ifelse(Sample_ID == "36-873_08092021", "August", "September"),
         MicrocystisBiomassFrac_ugCL = MicrocystisBiomassFrac*1000*0.22)

blah <- MicrocystisFracs %>% 
  group_by(Sample_ID) %>% 
  summarize(tot = sum(value))

Microbind <- rbind(L200, L100, L20)

pdf(file.path(here::here("03_Figures"),
              "16_FigS7_MicrocystisSizeFractions.pdf"), width = 8, height = 6)
ggplot(MicrocystisFracs, aes(y = MicrocystisBiomassFrac_ugCL, x = Month, fill = name)) +
  geom_bar(position = "stack", stat = "identity", color = "black")+
  ylab(expression(paste(italic("Microcystis "), "sp. biomass (µg C ",L^-1,")")))+
  guides(fill = guide_legend(title = "Size fraction")) +
  scale_fill_manual(values = c("#238b45", "#66c2a4", "#b2e2e2","#edf8fb" )) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16, color = "black"),
        text = element_text(family = "serif"))
dev.off()

# % below 60
# 80%
dim(dat[dat$Dcol < 60 & dat$Sample_ID == "36-873_08092021",])[1]/
  dim(dat[dat$Sample_ID == "36-873_08092021",])[1]

# 58%
dim(dat[dat$Dcol < 60 & dat$Sample_ID == "36-873_09102021",])[1]/
  dim(dat[dat$Sample_ID == "36-873_09102021",])[1]