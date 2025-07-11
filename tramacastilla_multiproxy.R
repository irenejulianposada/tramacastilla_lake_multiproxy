packages <- c("readxl", "dplyr", "tidyverse", "patchwork", "rioja", "vegan", "stringr",
              "janitor", "Bchron", "RRatepol", "scales", "tidyr", "lattice", "ggplot2",
              "usethis", "cowplot", "tapas", "nlme", "paletteer", "mgcv", "factoextra", "ggrepel",
              "grid", "corrplot", "psych", "tidypaleo", "riojaPlot", "scales", "ggpattern") #RRatepol only works with R 4.2.0
lapply(packages, library, character.only = TRUE)
usethis::use_git()

###### 0. Load data ####
load("E:/Saco/IJP/2_code/tramacastilla_chapter2/tram20_multiproxy.RData")

###### 1. Load ages and depths ####
#load("E:/Saco/IJP/2_code/tramacastilla_chrono/tram20_depth_age.RData")
tram20_depth_ages <- read.table("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM20/Tram20_composite_depth_ages_v2.csv", header = TRUE, sep = ",")
tram20_depth_ages_sed_rate <- tram20_depth_ages %>%
  rename("age" = "median") %>%
  select(-depth) %>%
  rename("depth" = "composite_depth") %>%
  mutate(depth_diff = c(NA, diff(depth)),
         age_diff = c(NA, diff(age)),
         sed_rate = age_diff/depth_diff)

###### 2. Tramacastilla charcoal ######
##### 2.1. Data upload ----
setwd("E:/Saco/IJP/1_data/Tramacastilla/charcoal/input data")
charcoal.tram <-read.table("Tramacastilla.csv", h=T, dec=".", sep=";")
samples_names <- read.table("Samples_names.csv", h=T, sep = ";")
Tramacastilla_age <- read.table("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM20/TRAM20_composite_depth_ages_v2.csv", h=T, dec=".", sep=",") %>%
  select(SampleId, median, composite_depth) %>%
  rename("age" = "median",
         "depth" = "composite_depth") %>%
  mutate(age = as.numeric(age),
         depth = as.numeric(depth)) %>%
  mutate(SampleId = str_replace_all(SampleId, "TRAM20-1B-", "TRAM-")) %>%
  filter(SampleId %in% samples_names$SampleId) %>% #para quedarme con secciones 3U y 4U
  mutate_all(~replace(., is.na(.), 0))
length(unique(charcoal.tram$SampleId)) #number of samples included
#I need to add samples that are not here because i should have 202 and i have 194
charcoal.tram.202 <- full_join(charcoal.tram, samples_names, by = "SampleId") %>%
  mutate(ImageType = ifelse(is.na(ImageType), "GLOBAL", ImageType)) #i change ImageType to say GLOBAL for these new samples i added
length(unique(charcoal.tram.202$SampleId)) #make sure I have all 202 samples

#now i want to duplicate this new rows so i can have one row for GLOBAL and other for SEEDLE
samples_no_data <- charcoal.tram.202 %>%
  filter(is.na(ImageId)) %>%
  mutate(ImageType = "SEEDLE")

charcoal.tram.all <- rbind(charcoal.tram.202, samples_no_data) %>%
  mutate_all(~replace(., is.na(.), 0))
length(unique(charcoal.tram.all$SampleId)) #make sure I have all 202 samples

Gchar <- subset(charcoal.tram.all, ImageType=="GLOBAL")
Schar <- charcoal.tram.all %>%
  subset(ImageType=="SEEDLE") %>%
  mutate(NPath.WLRatio = ifelse(NPath.WLRatio > 1, 1/NPath.WLRatio, NPath.WLRatio)) %>% #change values of WLRatio>1 and do 1/that value due to software mistake
  mutate(NofSeedles = as.numeric(NofSeedles)) %>%
  mutate(DPI.ProjAreaObj = as.numeric(DPI.ProjAreaObj)) %>%
  mutate(ObjectId = as.numeric(ObjectId))

length(unique(Gchar$SampleId)) #make sure I have 202 samples for GLOBAL
length(unique(Schar$SampleId))#make sure I have 202 samples for SEEDLE

##### 2.2. I want to make a new data frame only with charcoal particles larger than 22500 ----
#First i want to see the frequencies of area
hist(Schar$DPI.ProjAreaObj, breaks = seq(0, 850000, 5000))
summary(Schar$DPI.ProjAreaObj)
# 25% of the particles are below 13496 um2, which is a sieve diameter of 116 um, and mean 35713
# we are going to filter particles so we remove all that are below 110 um, that we assume is the noisy part of the signal
#I want to remove those charcoal with an area <110 um
seedles <- charcoal.tram.all %>%
  filter(ImageType == "SEEDLE" & DPI.ProjAreaObj > 12100) %>%
  mutate(NPath.WLRatio = ifelse(NPath.WLRatio > 1, 1/NPath.WLRatio, NPath.WLRatio)) %>% #change values of WLRatio>1 and do 1/that value due to software mistake
  as.data.frame() %>%
  mutate(DPI.ProjAreaObj = as.numeric(DPI.ProjAreaObj))

#other data frame with two columns, SampleId and CharCount that is the number of particles that meet the requirements
counting <- seedles %>%
  group_by(SampleId) %>%
  summarise(CharCount = n())

#other data frame with two columns, SampleId and sumCharArea that is the sum of area of particles that meet the requirements
totalarea <- seedles %>%
  group_by(SampleId) %>%
  summarise(sumCharArea = sum(DPI.ProjAreaObj, na.rm = TRUE))

summary(totalarea$sumCharArea)
hist(totalarea$sumCharArea, breaks = seq(0, 4500000, 50000))

seedles.age <- seedles %>%
  merge(Tramacastilla_age, by = "SampleId") %>%
  group_by(SampleId) %>%
  mutate(Color = ifelse(median(NPath.WLRatio) > 0.5, "green", "yellow")) %>%
  distinct() %>%
  ungroup() %>%
  as.data.frame()

#### 2.2.1. Create a data frame for plotting ----
# Create a file with one value per sample: count=length sampleID, W:L= mean W:L, area= sum ProjArea
CharTotS_G.110 <- merge(counting, totalarea)
CharTotS_G.110.age <- merge(CharTotS_G.110, Tramacastilla_age, by = "SampleId")

#to check if area and charcoal number are correalted
cor.test(CharTotS_G.110.age$sumCharArea,CharTotS_G.110.age$CharCount, method="spearman", exact = FALSE) #87% of the data is correlated
#to check if that correlation is linnear or not
summary(lm(CharTotS_G.110.age$sumCharArea~CharTotS_G.110.age$CharCount)) #78% is the amount of data linnearly correlated

#i join these two data frames to have all good charcoal particles with CharCount column
CharTotS_S.110 <- seedles %>%
  group_by(SampleId) %>%
  summarise(meanCharWL = mean(NPath.WLRatio, na.rm = TRUE),
            medianCharWL = median(NPath.WLRatio, na.rm = TRUE))

CharTotS.110 <- CharTotS_G.110 %>%
  full_join(CharTotS_S.110, by = "SampleId") %>%
  full_join(samples_names, by = "SampleId") %>% #i add samples that are missing
  mutate_all(~replace(., is.na(.), 0))

#if i want to order the samples i have to add ages
CharTotS.110.age <- full_join(CharTotS.110, Tramacastilla_age, by = "SampleId")
#write.table(CharTotS.110.age, "E:/Saco/IJP/1_data/Tramacastilla/charcoal/Charcoal metric Tramacastilla 3U-4U.csv", dec=".", sep=",", quote=F, row.names=F)

#### 2.2.2. Graphic representation ----
seedles %>% #min, max, mean of WLRatio
  ggplot(aes(x= 1, y = NPath.WLRatio)) +
  geom_boxplot(outlier.shape = 20, width = 1) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  theme_classic()

#### 2.2.3. CharCount & Area; W/L Ratio ----
CharTotS.110.age %>%
  ggplot(aes(x = age, y = CharCount)) +
  geom_line() +
  xlab("Age (cal yr BP)") +
  ylab("Number of charcoals") +
  ggtitle("Burnt biomass") +
  theme_classic() +
  scale_x_reverse(breaks = seq(500, 12000, 500))

#boxplot W/L Ratio with age and diferent colours according to median W/L ratio value
seedles.age %>%
  ggplot(aes(group = age, x = age, y = NPath.WLRatio)) +
  geom_boxplot(aes(colour = Color), show.legend = FALSE, width = 4) +
  labs(x= "Age (cal yr BP)", y = "W/L Ratio", title = "Type of fuel") +
  theme_classic() +
  scale_x_reverse(limits = c(max(seedles.age$age), min(seedles.age$age)), breaks = seq(500, 12000, by = 500)) +
  scale_colour_manual(values = c("yellow" = "#FFD700", "green" = "#548B54")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "red") 
#+ ggsave("E:/Saco/IJP/1_data/Tramacastilla/charcoal/figures/Type of fuel_W.L.Ratio_boxplot.pdf", height = 8, width = 15) #height = 5, width = 9.375 para ver las letras mas grandes

#other way of plotting W/L Ratio
W.L.ratio <- seedles.age %>%
  group_by(age, depth) %>%
  summarise(Woody = sum(NPath.WLRatio > 0.5),
            Grasslands = sum(NPath.WLRatio < 0.5)) %>%
  pivot_longer(cols = c("Woody", "Grasslands"), names_to = "fuel", values_to = "count")

W.L.ratio %>%
  ggplot (aes(x = age/1000, y = count, fill = fuel, color = fuel)) +
  geom_area(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Age (ka BP)", y = "Number of particles", fill = "Type of fuel") +
  theme_classic() +
  scale_x_reverse(limits = c(max(W.L.ratio$age)/1000, min(W.L.ratio$age)/1000),
                  breaks = seq(0.5, 12, by = 0.5),
                  labels = ifelse(seq(0.5, 12, by = 0.5) %% 1 == 0, as.character(seq(0.5, 12, by = 0.5)), "")) +
  scale_fill_manual(values = c("Woody" = "darkgreen", "Grasslands" = "#FFD700")) +
  scale_color_manual(values = c("Woody" = "darkgreen", "Grasslands" = "#FFD700"))

#Charcoal number and area
CC.110 <- CharTotS.110.age %>%
  ggplot(aes(x = age, y = CharCount)) +
  geom_line() +
  labs(x = "Age (cal yr BP)", y = "Charcoal number") +
  scale_x_reverse(breaks = seq(500, 12000, 500)) +
  scale_y_continuous(breaks = seq(min(CharTotS.110.age$CharCount),max(CharTotS.110.age$CharCount), by = 20)) +
  theme_classic() +
  geom_smooth(method = "loess", se = FALSE)

SCA.110 <- CharTotS.110.age %>%
  ggplot(aes(x = age, y = sumCharArea)) +
  geom_line() +
  labs(x = "Age (cal yr BP)", y = "Charcoal area (µm2)") +
  scale_x_reverse(breaks = seq(500, 12000, 500)) +
  scale_y_continuous(breaks = seq(min(CharTotS.110.age$sumCharArea),max(CharTotS.110.age$sumCharArea), by = 500000)) +
  theme_classic() + 
  geom_smooth(method = "loess", se = FALSE)

plot_grid(CC.110,SCA.110,ncol = 1, nrow = 2, align="v") #+
#ggsave("E:/Saco/IJP/1_data/Tramacastilla/charcoal/figures/Charcoal number and area.pdf", height = 8, width = 15)

###### 3. CHAR (charcoal particles · cm-2 · yr-1) ####
##### 3.1 Input data ####
mass.density <- read.table("E:/Saco/IJP/1_data/Tramacastilla/charcoal/input data/Tramacastilla_charcoal_mass_density.csv", sep = ";", header = T) #dataframe with volume and mass information of charcoal samples
counting #dataframe with sampleid and number of charcoal particles found in each sample
CharTotS.110.age #dataframe with sampleid and all charcoal metrics

sed.rate.charcoal <- TRAM20_age_depth_sed_rate %>% #dataframe with sedimentation rate of charcoal samples
  mutate(SampleId = str_replace_all(SampleId, "TRAM20-1B-", "TRAM-")) %>% #change SampleId so it is the same than in CHAR 
  mutate (across(c(depth_diff, age_diff, sed_rate), ~ ifelse(row_number() == 1 & is.na(.x), .x[2], .x))) %>% #we assume that sedimentation rate of first sample is the same than the second
  filter(SampleId %in% mass.density$SampleId)

##### 3.2. Calculate CHAR ####
#First calculate CHAR based on the number of particles
CHAR_number <- counting %>%
  full_join(mass.density, by = "SampleId") %>%
  select(-3) %>%
  rename("Mass" = "Mass..g.",
         "Density" = "Density..dry.sample...g.cm3.") %>%
  mutate_at("CharCount", ~replace(., is.na(.), 0)) %>% #NAs in CharCount are samples with no charcoal, so 0
  left_join(sed.rate.charcoal %>%
              select(SampleId, age, sed_rate, composite_depth) %>%
              rename("depth" = "composite_depth"), by = "SampleId") %>%
  relocate("age", .after = "SampleId") %>%
  mutate(age = as.numeric(age)) %>%
  mutate("CHAR.tram" = (CharCount/Mass) * Density * sed_rate) %>%
  mutate_all(~replace(., is.na(.), 0))

#Calculate CHAR based on the total area of charcoal particles
CHAR <- CharTotS.110.age %>%
  select(SampleId, sumCharArea, age, depth) %>%
  full_join(mass.density, by = "SampleId") %>%
  select(-"Density..wet.sample...g.cm3.") %>%
  rename("Mass" = "Mass..g.",
         "Density" = "Density..dry.sample...g.cm3.") %>%
  left_join(select(sed.rate.charcoal, SampleId, sed_rate), by = "SampleId") %>%
  relocate("age", .after = "SampleId") %>%
  mutate("CHAR.tram" = ((sumCharArea/Mass) * Density * sed_rate)/10e6) #mm2·cm-2·yr-1

#Correlation between CHAR and sedimentation rate
cor(CHAR$sed_rate, CHAR$CHAR.tram, method = 'pearson') #cor = 0.36 -> moderate correlation

##### 3.3. Figure CHAR ####
pCHAR <- CHAR %>%
  ggplot(aes(x = age/1000, y = CHAR.tram)) +
  geom_col(width = 0.05) +
  ylab("CHAR (charcoal area mm2·cm-2·yr-1)") +
  xlab("Age (ka BP)") +
  scale_x_reverse(limits = c(max(CHAR$age)/1000, min(CHAR$age)/1000),
                  breaks = seq(0.5, 12, by = 0.5),
                  labels = ifelse(seq(0.5, 12, by = 0.5) %% 1 == 0, as.character(seq(0.5, 12, by = 0.5)), "")) +
  theme_classic() #+
#ggsave("E:/Saco/IJP/1_data/Tramacastilla/charcoal/figures/CHAR_area.png", height = 5, width = 8)

##### 3.4. Figure 6 Cap2 #####
line_depths <- c(449, 497, 537, 620, 719, 765, 782, 822, 847)
zone_labels <- paste0("Z-", 10:1)
#add lines for coniss zones
add_linesy <- function(plot) {plot +
    geom_hline(yintercept = line_depths, color = "black", linetype = "dashed", size = 0.3, alpha = 0.7)}

fig6.1 <- CHAR %>%
  ggplot(aes(x = CHAR.tram, y = depth)) +
  geom_col(width = 1.2, orientation = "y") +
  xlab("CHAR \n(mm2·cm-2·yr-1)") +
  ylab("Depth (cm)") +
  scale_y_reverse(limits = c(880, 380),
                  breaks = seq(380, 880, by = 20) ,
                  labels = ifelse(seq(380, 880, by = 20) %% 40 == 0, as.character(seq(380, 880, by = 20)), "")) +
  scale_x_continuous(expand = c(0.01, 0))+
  theme_classic() +
  annotate("text", y = 800, x = 0, label = "No data", hjust = -1, color = "black")

fig6.1_coniss <- fig6.1 %>%add_linesy()

fig6.2 <- W.L.ratio %>%
  rename("Type of fuel" = "fuel") %>%
  ggplot(aes(x = count, y = depth, fill = `Type of fuel`, color = `Type of fuel`)) +
  geom_area(alpha = 0.6, position = "stack", orientation = "y") +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, aes(color = `Type of fuel`), linewidth = 1, alpha = 0.8, orientation = "y") +
  scale_x_continuous(trans = "sqrt",  # transformation with squared root
                     breaks = c(0, 5, 10, 20, 30, 50),
                     expand = c(0.01, 0)) +
  labs(y = "", x = "Number of particles \n(sqrt scale)") +
  theme_classic() +
  scale_y_reverse(limits = c(880, 380),
                  breaks = seq(380, 880, by = 20) ,
                  labels = ifelse(seq(380, 880, by = 20) %% 40 == 0, as.character(seq(380, 880, by = 20)), "")) +
  scale_fill_manual(values = c("Woody" = "darkgreen", "Grasslands" = "#FFD700")) +
  scale_color_manual(values = c("Woody" = "darkgreen", "Grasslands" = "#FFD700")) +
  annotate("text", y = 450, x = 120, label = "peak = 59", vjust = -1, color = "red") +
  annotate("text", y = 820, x = 10, label = "No data", vjust = -1, color = "black")

fig6.2_coniss <- fig6.2 %>%
  add_linesy() +
  annotate("text", y = (c(380, line_depths) + c(line_depths, 880)) / 2,  # Puntos medios
           x = Inf, label = zone_labels, hjust = 1.5, vjust = 0.5, size = 3, color = "black")

fig6.1_coniss/fig6.2_coniss +
  plot_layout(heights = c(1, 1), ncol = 2, nrow = 1, guides = "collect") +
  theme(plot.margin = margin(0, 0, 0, 0, "pt")) #save pdf A5

###### 4. Tramacastilla pollen ####
tram20_pollen <- read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/pollen/Conteo_TRAM.xlsx", col_names = FALSE) %>% #Pollen counts
  t() %>%
  as.data.frame() %>%
  row_to_names(row_number = 1) %>%
  mutate(across(-c(1, 2), as.numeric)) %>%
  left_join(tram20_depth_ages, by = "SampleId") %>%
  select(-c("Section and cm", "Taxon", "min", "max", "mean", "depth")) %>%
  relocate(c("composite_depth", "median"), .after = SampleId) %>%
  rename("age" = "median",
         "depth" = "composite_depth") %>%
  select(where(~ !all(is.na(.)))) #exclude taxa that do not appear in these samples

summary(tram20_pollen$`Terrestrial pollen sum`)
sd(tram20_pollen$`Terrestrial pollen sum`)
summary(tram20_pollen$`No of taxa`)
sd(tram20_pollen$`No of taxa`)

#plot pollen type-richness
tram20_pollen %>%
  ggplot(aes(x = age/1000, y = `No of taxa`)) +
  geom_point() +
  scale_x_reverse(limits = c(25, 0.863),
                  breaks = seq(0.5, 25, by = 0.5),
                  labels = ifelse(seq(0.5, 25, by = 0.5) %% 1 == 0, as.character(seq(0.5, 25, by = 0.5)), "")) +
  theme_classic() +
  geom_smooth(method = "loess")

pollen_group <- read.table("E:/Saco/IJP/1_data/Tramacastilla/pollen/pollen_groups.csv", header =  TRUE, sep = ",")
length(unique(pollen_group$pollen_group))
length(unique(pollen_group$pollen_group_general))
length(unique(pollen_group$pollen_group_discussion))

#Long format with taxa and percentages
tram20_pollen_long <- tram20_pollen %>%
  select(-c(`Evergreen trees`:Trees, Shrubs, Herbs, Hygrophytes, Hydrophytes:Spores, `Charcoal <150`)) %>%#first i remove columns that correspond to taxa groups
  pivot_longer(cols = c("Abies":"Undetermined"),
               values_to = 'pollen_counts',
               names_to = 'Taxon') %>%
  mutate_at("pollen_counts", ~replace(., is.na(.), 0))

#Long format with pollen groups
tram20_pollen_groups_long <- tram20_pollen_long %>%
  left_join(pollen_group, by = "Taxon") 

#Long format with pollen groups 2 and percentages
tram20_pollen_groups_general_long_perc <- tram20_pollen_groups_long %>%
  group_by(SampleId, pollen_group_general) %>%
  summarise(total_pollen = sum(pollen_counts, na.rm = TRUE)) %>%
  left_join(tram20_pollen_groups_long %>%
              select(SampleId, depth, age, Lycopodium, `Terrestrial pollen sum`) %>%
              distinct(), by = "SampleId") %>%
  mutate(percentage = (total_pollen/`Terrestrial pollen sum`)*100) %>% #compute group percentages in relation to terrestrial pollen sum
  select(SampleId, depth, age, Lycopodium, pollen_group_general, total_pollen, `Terrestrial pollen sum`, percentage) %>%
  ungroup()

#wide format with pollen groups 2 and percentages
tram20_pollen_groups_general_perc <- tram20_pollen_groups_general_long_perc %>%
  select(age, depth, pollen_group_general, percentage) %>%
  pivot_wider(names_from = pollen_group_general, values_from = percentage, values_fill = 0) %>%
  select(age, depth, Abies, Pinus, Juniperus, Fagus, Betula, Corylus, `Deciduous Quercus`, Tilia, Ulmus, `Evergreen trees`,
         `Other deciduous trees`, Ephedra, Shrubs, Poaceae, `Cerealia type`, Apiaceae, Artemisia, Cichorioideae, Caryophyllaceae,
         Chenopodiaceae, Helianthemum, Plantago, Ranunculaceae, Thalictrum, Filipendula, Potentilla, Sanguisorba, Rumex,
         Scrophulariaceae, Urticaceae, Urtica, `Other herbs`, Hygrophytes, Hydrophytes, Undetermined,
         everything()) %>% #keep order
  arrange(age) %>%
  as.data.frame()

#taxa percentages
tram20_pollen_long_perc <- tram20_pollen_long %>%
  group_by(SampleId, Taxon) %>%
  mutate(percentage = (pollen_counts/`Terrestrial pollen sum`)*100) %>%
  ungroup()

###### 5. PAR ######
tram20_density_pollen <- read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/density/TRAM20/Densidad TRAM20-1B-3U, 4U, 5U.xlsx", col_names = T) %>%
  mutate(SampleId = paste0(SONDEO, "-", cm)) %>%
  filter(SampleId %in% tram20_pollen$SampleId) %>%
  rename("Density" = "DENSIDAD MUESTRA SECA (g/cm3)") %>%
  select(SampleId, Density)

pollen_samples_mass <- read.table("E:/Saco/IJP/1_data/Tramacastilla/pollen/inventario_polen_lab_TRAMACASTILLA.csv",  header = TRUE, sep = ",") %>%
  select(MUESTRA, PESO..g.) %>%
  mutate(across(everything(), ~ str_replace(., "^PYC-", ""))) %>%
  mutate(across(everything(), ~ str_replace(., ", ", "-"))) %>%
  mutate(across(everything(), ~ str_replace(., " cm", ""))) %>%
  rename("SampleId" = "MUESTRA",
         "Mass" = "PESO..g.") %>%
  mutate(Mass = as.numeric(Mass))

sed.rate.pollen <- TRAM20_age_depth_sed_rate %>% #dataframe with sedimentation rate of charcoal samples
  mutate (across(c(depth_diff, age_diff, sed_rate), ~ ifelse(row_number() == 1 & is.na(.x), .x[2], .x))) %>% #we assume that sedimentation rate of first sample is the same than the second
  filter(SampleId %in% tram20_pollen$SampleId)

tram20_PAR <- tram20_pollen %>%
  select(-c(`Evergreen trees`:Trees, Shrubs, Herbs, Hygrophytes, Hydrophytes:Spores, `Charcoal <150`)) %>%
  left_join(tram20_density_pollen, by = "SampleId") %>%
  left_join(pollen_samples_mass, by = "SampleId") %>%
  left_join(select(sed.rate.pollen, SampleId, sed_rate), by = "SampleId") %>%
  mutate(age = as.numeric(age)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(across(!c(SampleId, depth, age, `Lycopodium`, Density, Mass, sed_rate),
                ~ (. / ((`Lycopodium` * Mass)/41696)) * Density * sed_rate,
                .names = "PAR_{.col}")) %>% #41696 lycopodium spores  in each sample
  select(c(depth, age, PAR_Abies:`PAR_Undetermined`)) %>%
  rename_with(~ gsub("^PAR_", "", .x)) %>% #remove prefix
  select(where(~ sum(., na.rm = TRUE) > 0)) #only taxons that appear in these samples

#sum of terrestrial PAR in each sample
#see relationship between PAR (terrestrial vegetation) and CHAR
PAR_sum <- tram20_PAR %>% mutate(sum = rowSums(across(Abies:Viola))) %>%
  left_join(CHAR %>% select(age, CHAR.tram), by = "age")

summary(PAR_sum$sum)
sd(PAR_sum$sum)

PAR_sum %>% ggplot(aes(x = age, y=sum)) +
  geom_col() +
  scale_x_reverse() +
  theme_classic()

PAR.CHAR <- PAR_sum %>%
  na.omit()

#to check if PAR and CHAR are correalted
cor.test(PAR.CHAR$sum, PAR.CHAR$CHAR.tram, method="spearman", exact = FALSE) #61% of the data is correlated
#to check if that correlation is linnear or not
summary(lm(PAR.CHAR$sum~PAR.CHAR$CHAR.tram)) #9% is the amount of data linnearly correlated

PAR.CHAR %>%
  ggplot(aes(x = CHAR.tram, y = sum)) + 
  geom_point(shape = 21, fill = "gray50", color = "black", size = 4, alpha = 0.5) +
  geom_smooth(method = "lm", col = "red") +
  xlab("CHAR (particles/cm2 yr)") + 
  ylab("total PAR (pollen grains/cm2 yr)") +
  ggtitle("PAR vs. CHAR") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size  =  16))

##### 5.1. PAR grouped #####
tram20_PAR_groups <- tram20_PAR %>%
  pivot_longer(cols = Abies:Undetermined, names_to = "Taxon", values_to = "PAR_value") %>% #long format
  left_join(pollen_group, by = "Taxon") %>% #Include pollen groups
  group_by(depth, age, pollen_group_general) %>%
  summarise(PAR_sum = sum(PAR_value, na.rm = TRUE)) %>% #calculate PAR sum for each group
  pivot_wider(names_from = pollen_group_general, values_from = PAR_sum, values_fill = 0) %>% #back to wide format
  select(age, depth, Abies, Pinus, Juniperus, Fagus, Betula, Corylus, `Deciduous Quercus`, Tilia, Ulmus,
         `Evergreen trees`, `Other deciduous trees`, Ephedra, Shrubs, Poaceae, `Cerealia type`, Apiaceae, Artemisia,
         Cichorioideae, Caryophyllaceae, Chenopodiaceae, Helianthemum, Plantago, Ranunculaceae, Thalictrum,
         Filipendula, Potentilla, Sanguisorba, Rumex, Scrophulariaceae, Urtica, Urticaceae, `Other herbs`,
         Hygrophytes, Hydrophytes, Undetermined) %>% #keep order
  as.data.frame()

##### 5.1.1. Clustering #####
tram.dist.PARgroup <- vegdist(tram20_PAR_groups, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)
tram.chclust.PARgroup <- chclust(tram.dist.PARgroup, method="coniss")
bstick(tram.chclust.PARgroup)
#to see depths and ages of coniss zones
coniss.zones.10 <- cutree(tram.chclust.PARgroup, k = 10)
cluster_boundaries <- which(diff(coniss.zones.10)!=0)
coniss.depths <- tram20_PAR_groups$depth[cluster_boundaries] #depths = c(447, 495, 535, 618, 717, 763, 778, 818, 845)
coniss.ages <- tram20_PAR_groups$age[cluster_boundaries] #ages = c(2003, 3576, 4173, 6192, 11280, 12713, 14131, 14909, 21663)

##### 5.2. PAR grouped groups_discussion #####
tram20_PAR_groups_discussion <- tram20_PAR %>%
  pivot_longer(cols = Abies:Undetermined, names_to = "Taxon", values_to = "PAR_value") %>% #long format
  left_join(pollen_group, by = "Taxon") %>% #Include pollen groups
  group_by(depth, age, pollen_group_discussion) %>%
  summarise(PAR_sum = sum(PAR_value, na.rm = TRUE)) %>% #calculate PAR sum for each group
  pivot_wider(names_from = pollen_group_discussion, values_from = PAR_sum, values_fill = 0) %>% #back to wide format
  select(age, depth, Abies, Pinus, Juniperus, Fagus, Betula, Corylus, `Deciduous Quercus`, Tilia, Ulmus,
         `Evergreen trees`, `Other deciduous trees`, Shrubs, Poaceae, `Steppe taxa`, `Other herbs`, `Anthropogenic taxa`,
         Hygrophytes, Hydrophytes, Undetermined) %>% #keep order and exclude Aquatics, Ferns, Selaginella
  as.data.frame()

##### 5.3. PAR all pollen types and aquatics - Supp Figure 5 Cap 2#####
# Plot bar plot
zone_labels_pollen <- paste0("TRAMZ-", 10:1)
labels_points <- (c(380, line_depths) + c(line_depths, 880)) / 2 #for coniss labels
p.col.all.1 <- c(rep("forestgreen", times=18), rep("brown", times = 18), rep("orange", times=18))
pollen.plot.PAR1 <- strat.plot(tram20_PAR %>%
                                 select(Abies:`Caprifoliaceae t Centranthus`),
                               yvar=tram20_PAR %>%
                                 select(depth),
                               ylabel="Depth (cm)", y.tks=seq(380, 880, 20), y.rev=TRUE,
                               plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, exag=TRUE,
                               srt.xlabel=70, col.poly= p.col.all.1, col.poly.line = "black", col.bar = "black",
                               scale.percent=TRUE, x.pc.inc=200000, xSpace=0.009, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2) #ADD CLUSTER ZONES FROM PAR VALUES
addClustZone(pollen.plot.PAR1, tram.chclust.PARgroup, nZone=10, lwd=1.5, lty=2, col="grey25") #save pdf 17*11
text(x = par("usr")[1] + 1.15 * diff(par("usr")[1:2]),
     y = labels_points, labels = zone_labels_pollen, pos = 2, cex = 0.6, xpd = TRUE)

p.col.all.2 <- c(rep("orange", times=48), rep("#009292", times = 5), rep("#33E4FF", times = 5), "#838B8B")
pollen.plot.PAR2 <- strat.plot(tram20_PAR %>%
                                 select(Caryophyllaceae:Undetermined),
                               yvar = tram20_PAR %>%
                                 select(depth),
                               ylabel="Depth (cm)", y.tks=seq(380, 880, 20), y.rev=TRUE,
                               plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, exag = TRUE, srt.xlabel=70,
                               col.poly = p.col.all.2,  col.bar= "black",col.poly.line = "black", scale.percent=TRUE, x.pc.inc=200000,
                               xSpace=0.009, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2)
addClustZone(pollen.plot.PAR2, tram.chclust.PARgroup, nZone=10, lwd=1.5, lty=2, col="grey25") #save pdf 17*11
text(x = par("usr")[1] + 1.15 * diff(par("usr")[1:2]),
     y = labels_points, labels = zone_labels_pollen, pos = 2, cex = 0.6, xpd = TRUE)

###### 6. Pollen percentages with all taxons - Supp Figure 4 Cap 2 ######
tram20_pollen_perc <- tram20_pollen_long_perc %>%
  select(depth, age, Taxon, percentage) %>%
  pivot_wider(names_from = Taxon,
              values_from = percentage)

# Figure 4 Supp Mat Cap 2
pollen.plot.all.1 <- strat.plot(tram20_pollen_perc %>%
                              select(Abies:`Caprifoliaceae t Centranthus`),
                            yvar = tram20_pollen_perc %>%
                              select(depth),
                            y.tks=seq(380, 880, 20), ylabel = "Depth (cm)", y.rev=TRUE,
                            plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, exag=TRUE,
                            srt.xlabel=70, col.poly= p.col.all.1, col.poly.line = "black", col.bar = "black", scale.percent=TRUE, x.pc.inc=5,
                            xSpace=0.001, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2) #ADD CLUSTER ZONES FROM PAR VALUES
addClustZone(pollen.plot.all.1, tram.chclust.PARgroup, nZone=10, lwd=1.5, lty=2, col="grey25") #save pdf 17*11
text(x = par("usr")[1] + 1.15 * diff(par("usr")[1:2]),
     y = labels_points, labels = zone_labels_pollen, pos = 2, cex = 0.6, xpd = TRUE)

pollen.plot.all.2 <- strat.plot(tram20_pollen_perc %>%
                                  select(Caryophyllaceae:Undetermined),
                                yvar = tram20_pollen_perc %>%
                                  select(depth),
                                y.tks=seq(380, 880, 20), ylabel = "Depth (cm)", y.rev=TRUE,
                                plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, exag=TRUE,
                                srt.xlabel=70, col.poly= p.col.all.2, col.poly.line = "black", col.bar = "black", scale.percent=TRUE, x.pc.inc=5,
                                xSpace=0.001, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2) #ADD CLUSTER ZONES FROM PAR VALUES
addClustZone(pollen.plot.all.2, tram.chclust.PARgroup, nZone=10, lwd=1.5, lty=2, col="grey25") #save pdf 17*11
text(x = par("usr")[1] + 1.15 * diff(par("usr")[1:2]),
     y = labels_points, labels = zone_labels_pollen, pos = 2, cex = 0.6, xpd = TRUE)

###### 7. Pollen and CHAR ######
##### 7.1. Fig 4 Cap 2 - Pollen, CHAR and sedimentation rate #####
#First i want to generate a dataframe with CHAR and sedimentation rate
tram20_CHAR_pollen_groups_general_perc <- tram20_pollen_groups_general_perc %>%
  full_join(CHAR %>% select(age, CHAR.tram), by = "age") %>%
  #mutate(across(where(is.numeric), ~ replace_na(., 0))) %>% #replace NAs
  select(-depth) %>%
  arrange(age) %>%
  mutate(across(c(CHAR.tram), ~ . *10)) %>% #BE CAREFUL WITH THESE. IT IS ONLY TO BETTER REPRESENT THE DATA
  rename("CHAR" = "CHAR.tram") %>%
  left_join(tram20_depth_ages_sed_rate %>%
              select(age, depth, sed_rate), by = "age") %>%
  mutate(across(sed_rate, ~ . /4)) %>% #BE CAREFUL WITH THESE. IT IS ONLY TO BETTER REPRESENT THE DATA
  relocate(depth, .after = "age") %>%
  relocate(sed_rate, .after = "depth")

#styles
c.sed.rate.line.1 <- c("black", rep(NA, times = 36))
ex.1 <- c(FALSE, rep(TRUE, times = 35), FALSE)
c.pollen.bar.1 <- c(NA, rep("black", times = 35), NA)
c.pollen.poly.1 <- c(NA, rep("forestgreen", times=11), rep("brown", times = 2), rep("orange", times=19), "#009292", "#33E4FF", "#838B8B", NA)
c.char.bar.1 <- c(rep(NA, times = 36), "black")

# Figure 4 manuscript Cap2
figure4_Cap2 <- strat.plot(tram20_CHAR_pollen_groups_general_perc %>%
                             select(sed_rate:CHAR) %>%
                             rename("sedimentation \n rate (yrs/cm)" = "sed_rate"),
                           yvar=tram20_CHAR_pollen_groups_general_perc %>%
                             select(depth),
                           ylabel="Depth (cm)", y.tks= seq(380, 880, 20), y.rev=TRUE,
                           plot.line=TRUE, plot.poly=FALSE, plot.bar=FALSE,
                           col.line = c.sed.rate.line.1, scale.percent=TRUE, exag = FALSE,
                           srt.xlabel=70, cex.xlabel = 0.8, x.pc.inc=5, cex.axis=0.6, 
                           xSpace=0.005, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2)
figure4_Cap2 <- strat.plot(tram20_CHAR_pollen_groups_general_perc %>%
                             select(sed_rate:CHAR),
                           yvar=tram20_CHAR_pollen_groups_general_perc %>%
                             select(depth),
                           y.tks= FALSE, y.rev=TRUE,
                           plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, col.bar = c.pollen.bar.1, col.poly= c.pollen.poly.1,
                           col.poly.line = c.pollen.bar.1,
                           scale.percent=TRUE, exag = ex.1, x.names="", x.axis=FALSE, 
                           xSpace=0.005, x.pc.lab=FALSE, x.pc.omit0=TRUE, las=2, add = TRUE)
figure4_Cap2 <- strat.plot(tram20_CHAR_pollen_groups_general_perc %>%
                             select(sed_rate:CHAR),
                           yvar=tram20_CHAR_pollen_groups_general_perc %>%
                             select(depth),
                           y.tks= FALSE, y.rev=TRUE,
                           plot.line=FALSE, plot.poly=FALSE, plot.bar=TRUE, col.bar = c.char.bar.1, lwd.bar = 6,
                           scale.percent=TRUE, exag = ex.1, x.names="", x.axis=FALSE, x.pc.inc=5,
                           xSpace=0.005, x.pc.lab=FALSE, x.pc.omit0=TRUE, las=2, add = TRUE)
#abline(h = tram20_CHAR_pollen_groups_general_perc %>%               ## DOES NOT WORK WELL
#         filter(depth %in% coniss.depths) %>%
#         select(depth) %>%
#         as_vector(),
#       lty = 2, col = "grey25", lwd = 1.5)
text(x = par("usr")[1] + 1.05 * diff(par("usr")[1:2]),
     y = labels_points, labels = zone_labels_pollen, pos = 2, cex = 0.6, xpd = TRUE)
#save pdf A4 and divide x axis by CHAR/10 and multiply sed_rate by 4. CUIDADO QUE NO SALE BIEN EL CLUSTER

##### 7.2. Fig 5 Cap 2 - PAR, CHAR and sedimentation rate #####
tram20_CHAR_PAR_groups <- tram20_PAR_groups %>%
  full_join(CHAR %>% select(age, CHAR.tram), by = "age") %>%
  select(-depth) %>%
  arrange(age) %>%
  left_join(tram20_depth_ages_sed_rate %>%
              select(age, depth, sed_rate), by = "age") %>%
  relocate(depth, .after = "age") %>%
  relocate(sed_rate, .after = "depth") %>%
  mutate(across(CHAR.tram, ~ . *1000000)) %>% #BE CAREFUL WITH THESE. IT IS ONLY TO BETTER REPRESENT THE DATA
  mutate(across(sed_rate, ~ . *20000)) %>% #BE CAREFUL WITH THESE. IT IS ONLY TO BETTER REPRESENT THE DATA
  rename("CHAR" = "CHAR.tram")

figure5_Cap2 <- strat.plot(tram20_CHAR_PAR_groups %>%
                             select(sed_rate:CHAR) %>%
                           rename("sedimentation \n rate (yrs/cm)" = "sed_rate"),
                           yvar=tram20_CHAR_pollen_groups_general_perc %>%
                             select(depth),
                           ylabel="Depth (cm)", y.tks= seq(380, 880, 20), y.rev=TRUE,
                           plot.line=TRUE, plot.poly=FALSE, plot.bar=FALSE,
                           col.line = c.sed.rate.line.1, scale.percent=TRUE, exag = FALSE,
                           srt.xlabel=70, cex.xlabel = 0.8, x.pc.inc=200000, cex.axis=0.6, 
                           xSpace=0.008, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2)
figure5_Cap2 <- strat.plot(tram20_CHAR_PAR_groups %>%
                             select(sed_rate:CHAR),
                           yvar=tram20_CHAR_pollen_groups_general_perc %>%
                             select(depth),
                           y.tks= FALSE, y.rev=TRUE,
                           plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, col.bar = c.pollen.bar.1, col.poly = c.pollen.poly.1,
                           col.poly.line = c.pollen.bar.1, scale.percent=TRUE, exag = ex.1, x.names="", x.axis=FALSE,
                           srt.xlabel=70, x.pc.inc=200000, xSpace=0.008, x.pc.lab=FALSE, x.pc.omit0=TRUE, las=2, add = TRUE)
figure5_Cap2 <- strat.plot(tram20_CHAR_PAR_groups %>%
                             select(sed_rate:CHAR),
                           yvar=tram20_CHAR_pollen_groups_general_perc %>%
                             select(depth),
                           y.tks= FALSE, y.rev=TRUE,
                           plot.line=FALSE, plot.poly=FALSE, plot.bar=TRUE, col.bar = c.char.bar.1, lwd.bar = 6,
                           scale.percent=TRUE, exag = FALSE, x.names="", x.axis=FALSE,
                           srt.xlabel=70, x.pc.inc=200000, xSpace=0.008, x.pc.lab=FALSE, x.pc.omit0=TRUE, las=2, add = TRUE)
text(x = par("usr")[1] + 1.05 * diff(par("usr")[1:2]),
     y = labels_points, labels = zone_labels_pollen, pos = 2, cex = 0.6, xpd = TRUE)
#save pdf A4 and divide x axis of CHAR by 1000000 and sed_rate by 20000

###### 8. Temperature ~ CHAR ######
bins <- CHAR$age
Mendukilo <- read.csv2("E:/Saco/IJP/1_data/Mendukilo/Mendukilo.csv", header=TRUE, sep=",", stringsAsFactors = FALSE) %>%
  mutate_all(as.numeric) %>%
  mutate(`age..kyr.BP.` = `age..kyr.BP.`*1000) %>%
  rename("d13C.Mendukilo" = "d13C..permil.",
         "age" = "age..kyr.BP.")

Mendukilo_binning <- Mendukilo %>%
  mutate(binned_age = sapply(age, function(x) bins[which.min(abs(bins-x))]))

Mendukilo_binned <- Mendukilo_binning %>%
  group_by(binned_age) %>%
  summarise(d13C.Mendukilo = mean(d13C.Mendukilo, na.rm = TRUE)) %>%
  rename("age" = "binned_age") %>%
  ungroup() %>%
  as.data.frame()

temp.char <- Mendukilo_binned %>%
  left_join(CHAR %>%
              select(age, CHAR.tram), by = "age")
pairs(temp.char)
temp.char %>%
  ggplot(aes(x = CHAR.tram, y = d13C.Mendukilo)) + 
  geom_point(shape = 21, fill = "gray50", color = "black", size = 4, alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, col = "red") +
  xlab("CHAR (particles/cm2 yr)") + 
  ylab("d13C Mendukilo") +
  ggtitle("Temperature vs. CHAR") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size  =  16))

temp.char.gls <- gls(d13C.Mendukilo ~ CHAR.tram, data = temp.char)
summary(temp.char.gls)
cor(temp.char$d13C.Mendukilo, predict(temp.char.gls))^2 #R2 = 0.08

###### 9. XRF data ######
#import table with all composite depths, after removing cm that do not have good XRF values (they are in red in XRF sheet) 
tram20_depths_XRF <- read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM20/TRAM20_core_depth_v4_XRF.xlsx", col_names = TRUE) %>%
  select(c(2:8, 7)) %>%
  mutate(site = ifelse(site == "TRAMACASTILLA", "TRAM20")) %>%
  mutate(SampleId = paste(site, core, section, `section depth cm`, sep = "-")) %>% #concateno varias columnas para crear la columna SampleId
  relocate(SampleId, .before = site) %>%
  rename("composite_depth" = "composite depth (cm)") %>%
  select(SampleId, composite_depth)

tram20_LECO_3U_4U_5U <- read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/LECO/20220829_LECO_TRAM20-1B.xlsx", col_names = FALSE) %>%
  setNames(as.character(slice(., 11))) %>%
  slice(-c(1:11)) %>%
  separate(Id, into = c("rest", "Id"), sep = "PYC-") %>%
  separate(cm, into = c("cm", "other"), sep = "-") %>%
  mutate(SampleId = paste(Id, cm, sep = "-")) %>%
  mutate(across(`CT (%)`:`TOC (%)`, as.numeric)) %>%
  mutate("TC" = `CT (%)`,
         "TN" = `NT (%)`,
         "TS" = `ST (%)`,
         "TIC" = `TIC (%)`,
         "TOC" = `TOC (%)`,
         "TOC/TN" = `TOC (%)`/`NT (%)`) %>%
  select(SampleId, TC, TN, TS, TIC, TOC, `TOC/TN`) %>%
  as.data.frame() %>%
  mutate(across(TC:TOC, as.numeric))

tram20_LECO_1U_2U <- read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/LECO/LECO_TRAM20-1B_81221_81341_CNS_TIC_INFORME.xlsx", col_names = FALSE) %>%
  setNames(as.character(slice(., 11))) %>%
  slice(-c(1:11)) %>%
  separate(cm, into = c("cm", "other"), sep = "-") %>%
  mutate(SampleId = paste(Id, cm, sep = "-")) %>%
  mutate(across(`TC (%)`:`TS (%)`, as.numeric)) %>%
  mutate("TC" = `TC (%)`,
         "TN" = `TN (%)`,
         "TS" = `TS (%)`,
         "TIC" = `TIC (%)`,
         "TOC" = `TOC (%)`,
         "TOC/TN" = `TOC (%)`/`TN (%)`) %>%
  select(SampleId, TC, TN, TS, TIC, TOC,`TOC/TN`) %>%
  as.data.frame()

tram20_LECO <- bind_rows(tram20_LECO_1U_2U, tram20_LECO_3U_4U_5U)

tram20_MS_IPE <- read_xls("E:/Saco/IJP/1_data/Tramacastilla/MS_IPE/TRAM20/TRAM20-1B UWITEC.xls", col_names = FALSE) %>%
  setNames(as.character(slice(., 9))) %>%
  mutate(`SECT DEPTH` = as.numeric(`SECT DEPTH`)) %>%
  slice(-c(1:11)) %>%
  mutate(SampleId = paste0("TRAM20-1B-", `SECT NUM`, "-", trunc(`SECT DEPTH`))) %>%
  mutate("MS" = as.numeric(MS1)) %>%
  select(SampleId, MS)

tram20_XRF_LECO_MS <- read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/XRF_MS_BCN/TRAM20/XRF/TRAM20_1B_ALL_modified_preliminary_ages.xlsx", sheet = 1, col_names = TRUE) %>%
  select(-c(`depth composite (no sed loss, no compress., april 2021)`, `ageBP (CLAM model april 2021)`)) %>% #remove old depths
  separate(Spectrum, into = c("section", "rest"), sep = " ", extra = "merge") %>% #separate sample names and cm
  mutate(SampleId = paste(section, `depth cm`, sep = "-")) %>% #to make it coincide with tram20_depths_XRF
  full_join(read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/XRF_MS_BCN/TRAM20/XRF/TRAM20_1B_ALL_modified_preliminary_ages.xlsx", sheet = 2, col_names = TRUE) %>%
              separate(Spectrum, into = c("section", "rest"), sep = " ", extra = "merge") %>% #separate sample names and cm
              mutate(SampleId = paste(section, Sample/10, sep = "-")),
            by = "SampleId") %>%
  right_join(tram20_depths_XRF, by = "SampleId") %>%
  mutate("K" = `K -Ka Area`,
         "Si" = `Si-Ka Area`,
         "Ca/Ti" = ifelse(`Ca-Ka Area`/`Ti-Ka Area` > 1.5, 1.5, `Ca-Ka Area`/`Ti-Ka Area`), #to facilitate visualisation peaks in Ca/Ti with a value >1.5 will be reduced to 1.5
         "S" = `S -Ka Area`,
         "Br" = `Br-Ka Area`,
         "Fe/Mn" = `Fe-Ka Area`/`Mn-Ka Area`,
         "Rb/Zr" = `Rb-Ka Area`/`Zr-Ka Area`) %>%
  full_join(tram20_LECO, by = "SampleId") %>%
  full_join(tram20_MS_IPE, by = "SampleId") %>%
  select(SampleId, composite_depth, TN, TS, TIC, TOC, `TOC/TN`, MS, K, Si, `Ca/Ti`, S, Br, `Fe/Mn`, `Rb/Zr`)

tram20_XRF_LECO_MS %>%
  pivot_longer(cols = c("TN":"Rb/Zr"),
               values_to = 'values',
               names_to = 'elements') %>%
  mutate(facet_group = case_when(elements == "TN" ~ "TN",
                                 elements == "TS" ~ "TS",
                                 elements == "TIC" ~ "TIC",
                                 elements == "TOC" ~ "TOC",
                                 elements == "TOC/TN" ~ "TOC/TN",
                                 elements == "MS" ~ "MS",
                                 elements %in% c("K", "Si") ~ "K & Si",       # i want to have K and Si together in the grid
                                 elements == "S" ~ "S",
                                 elements == "Br" ~ "Br",
                                 elements == "Ca/Ti" ~ "Ca/Ti",
                                 elements == "Fe/Mn" ~ "Fe/Mn",
                                 elements == "Rb/Zr" ~ "Rb/Zr"),
         facet_group = fct_relevel (facet_group, "TN", "TS", "TIC", "TOC", "TOC/TN", "MS", "K & Si", "Ca/Ti", "S", "Br", "Fe/Mn", "Rb/Zr")) %>% 
  ggplot(aes(x = values, y = composite_depth, fill = elements)) +
  geom_point(aes(color = elements), size = 0.5) + geom_lineh(aes(color = elements)) +
  scale_y_reverse(breaks = seq(0, 900, 100)) +
  scale_x_continuous() +
  facet_geochem_gridh(vars(facet_group)) +
  labs(y = "Depth (cm)",
       x = "") +
  theme_paleo() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("darkblue", "brown", "red", "#6551CC", "#666666", "darkgrey",
                               "black", "#F76D5E", "purple", "#FF7F00", "black", "#EE4C97", "#F100F1"))
  #+ ggsave("E:/Saco/IJP/5_manuscripts/Cap 2-QSR/figures/LECO_MS_XRF.pdf", width = 8, height = 7)

##### 9.1. Figure 3 Cap 2 - XRF data with MS, LECO, XRF #####
#ordered: MS, K & Si, Ca/Ti, S, TIC, TOC, Inc/Coh, Br, Fe/Mn, Rb/Zr
tram20_XRF_LECO_MS %>%
  select(SampleId, composite_depth, MS) %>%
  pivot_longer(cols = "MS",
               values_to = 'values',
               names_to = 'elements') %>%
  na.omit() %>%
  ggplot(aes(x = values, y = composite_depth, fill = elements)) +
  geom_lineh(aes(color = elements)) +
  scale_y_reverse(limits = c(900, 0), breaks = seq(0, 900, 100), labels = seq(0, 900, 100)) +
  scale_x_continuous(n.breaks = 3) +
  facet_geochem_gridh(vars(elements)) +
  theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_color_manual(values = "#666666") +
  tram20_XRF_LECO_MS %>%
  select(c(SampleId, composite_depth, K, Si, S,`Ca/Ti`)) %>%
  pivot_longer(cols = c("K":"Ca/Ti"),
               values_to = 'values',
               names_to = 'elements') %>%
  mutate(facet_group = case_when(elements %in% c("K", "Si") ~ "K & Si",       # i want to have K and Si together in the grid
                                 elements == "S" ~ "S",
                                 elements == "Ca/Ti" ~ "Ca/Ti"),
         facet_group = fct_relevel (facet_group, "K & Si", "Ca/Ti", "S")) %>%
  ggplot(aes(x = values, y = composite_depth, fill = elements)) +
  geom_lineh(aes(color = elements)) +
  scale_y_reverse(breaks = seq(0, 900, 100)) +
  scale_x_continuous(n.breaks = 3) +
  facet_geochem_gridh(vars(facet_group)) +
  theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_color_manual(values = c("brown", "#6551CC", "purple", "#FF7F00")) +
  tram20_XRF_LECO_MS %>%
  select(SampleId, composite_depth, TIC, TOC) %>%
  pivot_longer(cols = c("TIC":"TOC"),
               values_to = 'values',
               names_to = 'elements') %>%
  mutate(elements = fct_relevel (elements, "TIC", "TOC")) %>%
  na.omit() %>%
  ggplot(aes(x = values, y = composite_depth, fill = elements)) +
  geom_lineh(aes(color = elements)) +
  scale_y_reverse(limits = c(900, 0), breaks = seq(0, 900, 100), labels = seq(0, 900, 100)) +
  scale_x_continuous(n.breaks = 3) +
  facet_geochem_gridh(vars(elements)) +
  theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_color_manual(values = c("#EE4C97", "#F100F1")) +
  tram20_XRF_LECO_MS %>%
  select(c(SampleId, composite_depth, Br, `Fe/Mn`, `Rb/Zr`)) %>%
  pivot_longer(cols = c("Br":"Rb/Zr"),
               values_to = 'values',
               names_to = 'elements') %>%
  mutate(facet_group = case_when(elements == "Br" ~ "Br",
                                 elements == "Fe/Mn" ~ "Fe/Mn",
                                 elements == "Rb/Zr" ~ "Rb/Zr"),
         facet_group = fct_relevel (facet_group, "Br", "Fe/Mn", "Rb/Zr")) %>%
  ggplot(aes(x = values, y = composite_depth, fill = elements)) +
  geom_lineh(aes(color = elements)) +
  scale_y_reverse(breaks = seq(0, 900, 100)) +
  scale_x_continuous(n.breaks = 3) +
  facet_geochem_gridh(vars(facet_group)) +
  theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_color_manual(values = c("darkblue", "red", "#FF7F00")) +
  plot_layout(widths = c(1, 3, 2, 4)) #+
  #ggsave("E:/Saco/IJP/5_manuscripts/Cap 2-QSR/figures/All files/LECO_MS_XRF_v5.pdf", width = 8, height = 7)

##### 9.2. Supplementary Figure 2 Cap 2 - all XRF, MS, LECO data #####
#Supp. Figure 2 - all XRF, LECO, MS
tram20_XRF_LECO_MS_SuppFig <- read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/XRF_MS_BCN/TRAM20/XRF/TRAM20_1B_ALL_modified_preliminary_ages.xlsx", sheet = 1, col_names = TRUE) %>%
  select(-c(`depth composite (no sed loss, no compress., april 2021)`, `ageBP (CLAM model april 2021)`)) %>% #remove old depths
  separate(Spectrum, into = c("section", "rest"), sep = " ", extra = "merge") %>% #separate sample names and cm
  mutate(SampleId = paste(section, `depth cm`, sep = "-")) %>% #to make it coincide with tram20_depths_XRF
  full_join(read_xlsx("E:/Saco/IJP/1_data/Tramacastilla/XRF_MS_BCN/TRAM20/XRF/TRAM20_1B_ALL_modified_preliminary_ages.xlsx", sheet = 2, col_names = TRUE) %>%
              separate(Spectrum, into = c("section", "rest"), sep = " ", extra = "merge") %>% #separate sample names and cm
              mutate(SampleId = paste(section, Sample/10, sep = "-")),
            by = "SampleId") %>%
  right_join(tram20_depths_XRF, by = "SampleId") %>%
  mutate("K" = `K -Ka Area`,
         "Si" = `Si-Ka Area`,
         "Ca" = `Ca-Ka Area`,
         "Ti" = `Ti-Ka Area`,
         "S" = `S -Ka Area`,
         "Br" = `Br-Ka Area`,
         "Fe" = `Fe-Ka Area`,
         "Mn" = `Mn-Ka Area`,
         "Rb" = `Rb-Ka Area`,
         "Zr" = `Zr-Ka Area`,
         "Pb" = `Pb-La Area`,
         "Sr" = `Sr-Ka Area`) %>%
  full_join(tram20_LECO, by = "SampleId") %>%
  full_join(tram20_MS_IPE, by = "SampleId") %>%
  select(SampleId, composite_depth, TN, TS, TIC, TOC, MS,
         K, Si, Ca, Ti, S, Br, Fe, Mn, Rb, Zr, Pb, Sr)

#plot
tram20_XRF_LECO_MS_SuppFig %>%
  select(SampleId, composite_depth, MS) %>%
  pivot_longer(cols = "MS",
               values_to = 'values',
               names_to = 'elements') %>%
  na.omit() %>%
  ggplot(aes(x = values, y = composite_depth, fill = elements)) +
  geom_lineh(aes(color = elements)) +
  scale_y_reverse(limits = c(900, 0), breaks = seq(0, 900, 100), labels = seq(0, 900, 100)) +
  scale_x_continuous(n.breaks = 3) +
  facet_geochem_gridh(vars(elements)) +
  theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_color_manual(values = "#666666") +
  tram20_XRF_LECO_MS_SuppFig %>%
  select(SampleId, composite_depth, TN, TS, TIC, TOC) %>%
  pivot_longer(cols = c("TN":"TOC"),
               values_to = 'values',
               names_to = 'elements') %>%
  mutate(elements = fct_relevel (elements, "TN", "TS", "TIC", "TOC")) %>%
  na.omit() %>%
  ggplot(aes(x = values, y = composite_depth, fill = elements)) +
  geom_lineh(aes(color = elements)) +
  scale_y_reverse(limits = c(900, 0), breaks = seq(0, 900, 100), labels = seq(0, 900, 100)) +
  scale_x_continuous(n.breaks = 3) +
  facet_geochem_gridh(vars(elements)) +
  theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_color_manual(values = c("#FF4040", "#53868B", "#EE4C97", "#F100F1")) +
  tram20_XRF_LECO_MS_SuppFig %>%
  select(c(SampleId, composite_depth, K, Si, Ca, Ti, S, Br, Fe, Mn, Rb, Zr, Pb, Sr)) %>%
  pivot_longer(cols = c("K":"Sr"),
               values_to = 'values',
               names_to = 'elements') %>%
  mutate(elements = fct_relevel (elements, "K", "Si", "S", "Ca", "Ti", "Br", "Fe", "Mn", "Rb", "Zr", "Pb", "Sr")) %>%
  ggplot(aes(x = values, y = composite_depth, fill = elements)) +
  geom_lineh(aes(color = elements)) +
  scale_y_reverse(breaks = seq(0, 900, 100)) +
  scale_x_continuous(n.breaks = 3) +
  facet_geochem_gridh(vars(elements)) +
  theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_color_manual(values = c("brown", "#6551CC", "purple", "#FF7F00", "#458B00", "darkblue",
                                "red", "#FF7F00", "#9932CC", "#FFB90F", "#CD1076", "#8B814C")) +
  plot_layout(widths = c(1, 4, 13)) #+
#ggsave("E:/Saco/IJP/5_manuscripts/Cap 2-QSR/supplementary material/Supp.Fig.2_MS_LECO_XRF.pdf", width = 10, height = 7)

##### 10. Diatoms #####
tram20_diatoms <- read.table("E:/Saco/IJP/1_data/Tramacastilla/diatoms/PYCTRAM_diatoms_v2.csv", header = TRUE, sep = ";") %>%
  select(-depth) %>%
  mutate(sample = str_replace_all(sample, "_", "-"),
         sample = str_replace_all(sample, "cm", ""),
         sample = paste0("TRAM20-1B-", sample)) %>%
  rename("SampleId" = "sample") %>%
  left_join(tram20_depth_ages %>%
              select(SampleId, composite_depth, median) %>%
              rename("age" = "median",
                     "depth" = "composite_depth"),
            by = "SampleId") %>%
  relocate("age", "depth", .after = "SampleId") %>%
  filter(age<15000) %>%
  select(where(~ !all(is.na(.)))) %>% #remove diatoms that do not appear in any sample
  mutate(across(everything(), ~replace_na(., 0))) %>%
  mutate(across(-SampleId, as.numeric)) %>%
  mutate(sum_diatoms = rowSums(across(Achnanthes.brevipes:Ulnaria.sp))) %>%
  mutate(richness = rowSums(select(., -c(SampleId, age, depth)) >0)) %>%
  relocate(sum_diatoms, .after = depth) %>%
  relocate(richness, .after = sum_diatoms)

#diatoms accumulation
summary(tram20_diatoms$sum_diatoms)
tram20_diatoms %>%
  ggplot(aes(x = age/1000, y = sum_diatoms)) +
  geom_point() +
  geom_line() +
  scale_x_reverse(breaks = seq(0, 15, 1)) +
  theme_classic()

#diatoms richness
summary(tram20_diatoms$richness)

# Fig. 7 - Selection of diatoms
inc_diatoms <- rep(5, ncol(tram20_diatoms %>%
                     select(Achnanthes.brevipes:Ulnaria.sp) %>%
                     select(where(~ sum(., na.rm = TRUE) > 5))))
ex.diat <- rep(TRUE, times = 55)
riojaPlot(tram20_diatoms %>%
            select(Achnanthes.brevipes:Ulnaria.sp) %>%
            select(where(~ sum(., na.rm = TRUE) > 5)), 
          tram20_diatoms %>%
            select(depth) %>%
            rename("Depth (cm)" = "depth"),
          yinterval = 20, ytks1=seq(380, 800, 20), scale.percent = TRUE, lwd.bar = 2, symb.cex = 8, symb.pch = 4,
          plot.exag = TRUE, srt.xlabel = 70, min.width.pc=0, x.pc.inc=inc_diatoms,
          cex.xaxis=0.7, cex.xlabel=0.8, las.xaxis=2)

diatoms_selection <- strat.plot(tram20_diatoms %>%
                                  select(Achnanthes.brevipes:Ulnaria.sp) %>%
                                  select(where(~ sum(., na.rm = TRUE) > 5)),
                                yvar=tram20_diatoms %>%
                                  select(depth),
                                ylabel="Depth (cm)", y.tks=seq(380, 800, 20), y.rev=TRUE, scale.percent = TRUE,
                                plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, col.bar="black", col.poly = "#1C86EE", col.poly.line = "black",
                                srt.xlabel=70, cex.xlabel = 1, cex.axis=0.8,
                                exag = ex.diat, x.pc.inc=5,
                                x.pc.omit0=TRUE, las=2, xSpace=0.007)
text(x = par("usr")[1] + 1.2 * diff(par("usr")[1:2]),
     y = labels_points, labels = zone_labels_pollen, pos = 2, cex = 0.6, xpd = TRUE)

#all diatoms for Supp Mat
inc_diatoms_all <- rep(5, ncol(tram20_diatoms %>%
                             select(Achnanthes.brevipes:Ulnaria.sp)))
ex.diatsupp <- rep(TRUE, times = 113)
riojaPlot(tram20_diatoms %>%
            select(Achnanthes.brevipes:Ulnaria.sp), 
          tram20_diatoms %>%
            rename("Depth (cm)" = "depth") %>%
            select(`Depth (cm)`),
          yinterval = 1000, ytks1=seq(380, 880, 20), scale.percent = TRUE, lwd.bar = 2, symb.cex = 8, symb.pch = 4,
          plot.exag = TRUE, srt.xlabel = 70, min.width.pc=0, x.pc.inc=5,
          cex.xaxis=0.7, cex.xlabel=0.7, las.xaxis=2, xSpace=0.002) #save pdf 17*11

strat.plot(tram20_diatoms %>%
             select(Achnanthes.brevipes:Ulnaria.sp),
           yvar=tram20_diatoms %>%
             select(depth),
           ylabel="Depth (cm)", y.tks=seq(380, 800, 20), y.rev=TRUE, scale.percent = TRUE,
           plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, col.bar="black", col.poly = "#1C86EE", col.poly.line = "black",
           srt.xlabel=70, cex.xlabel = 0.8, cex.axis=0.7,
           exag = ex.diatsupp, x.pc.inc=5,
           x.pc.omit0=TRUE, las=2, xSpace=0.005) #save pdf 17*11

###### 11. Discussion figure vertical #######
##### 11.1. Import sedaDNA data and create data frames for plotting #####
#### 11.1.1. Plant sedaDNA ####
plant_dna <- read.table("E:/Saco/IJP/1_data/Tramacastilla/sedaDNA_metabarcoding/TRAM21_plants_v3.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)

plant_rai <- plant_dna %>% #obtain dataframe
  mutate(across(.cols = starts_with("weightrep_sample.TRAM_1b_"),
                .fns = ~ if_any(starts_with("proprep_sample.TRAM_1b_"), ~ . <= 0.125) * 0 + 
                  if_all(starts_with("proprep_sample.TRAM_1b_"), ~ . > 0.125) * .)) %>% #i discard samples unreplicated
  mutate(across(starts_with("totread_sample."),
                ~. / sum(., na.rm = TRUE), .names = "norm_{col}")) %>% #normalize totread
  mutate(across(starts_with ("norm_totread_sample."),
                ~. * get(gsub("norm_totread_sample.", "weightrep_sample.", cur_column())),
                .names = "RAI_{sub('norm_totread_sample.', '', col)}")) %>% #totread*weightrep
  ungroup() %>%
  rename_with(~ gsub("RAI_.norm_totread_sample.", "RAI_", .), starts_with("RAI_.norm_totread")) %>% #rename
  rename_with(~ gsub("_rpt$", "", .), starts_with("RAI_")) %>%  #rename
  select(-c(starts_with("totread_sample"), starts_with("weightrep_sample"), starts_with("norm_totread_sample."))) %>%
  select(-starts_with("RAI_TRAM_C")) #remove controls

plant_rai_perc <- plant_rai %>%
  filter(!group %in% c("Bryophyte", "Hydrophyte", "Hygrophyte",
                       "Not native",  "Algae", "Positive control",
                       "Other plant", "Other tree")) %>% #remove some groups from percentages so that we obtain 
  mutate(across(starts_with("RAI_TRAM"), ~./sum(., na.rm = TRUE)*100))

plant_rai_long <- plant_rai %>%
  rename_with(~ gsub("RAI_", "", .), starts_with("RAI_")) %>%
  pivot_longer(cols = c("TRAM_1b_3U_11.12":"TRAM_1b_4U_109.110"),
               values_to = 'RAI',
               names_to = 'SampleId') %>%
  left_join(read.csv ("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM21/TRAM21_composite_depth_ages.csv") %>%
              rename("SampleId" = "SampleID",
                     "age" = "median") %>%
              select("SampleId", "age"),
            by = "SampleId")

plant_rai_perc_long <- plant_rai_perc %>%
  rename_with(~ gsub("RAI_", "", .), starts_with("RAI_")) %>%
  pivot_longer(cols = c("TRAM_1b_3U_11.12":"TRAM_1b_4U_109.110"),
               values_to = 'RAI',
               names_to = 'SampleId') %>%
  left_join(read.csv ("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM21/TRAM21_composite_depth_ages.csv") %>%
              rename("SampleId" = "SampleID",
                     "age" = "median") %>%
              select("SampleId", "age"),
            by = "SampleId")

#trees
trees_dna <- plant_rai %>%
  filter(scientific_name %in% c("Abies alba", "Betula", "Betulaceae", "Fagus",
                                "Fagus sylvatica", "Pinus", "Quercus 1")) %>%
  rename_with(~ gsub("RAI_", "", .), starts_with("RAI_")) %>%
  select(scientific_name, TRAM_1b_3U_11.12:TRAM_1b_4U_109.110) %>%
  mutate(scientific_name = case_when(scientific_name %in% c("Fagus", "Fagus sylvatica") ~ "Fagus sylvatica",
                           TRUE ~ scientific_name )) %>%
  mutate(group = paste0(scientific_name, "_dna")) %>%
  pivot_longer(cols = c("TRAM_1b_3U_11.12":"TRAM_1b_4U_109.110"),
               values_to = 'RAI',
               names_to = 'SampleId') %>%
  left_join(read.csv ("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM21/TRAM21_composite_depth_ages.csv") %>%
              rename("SampleId" = "SampleID",
                     "age" = "median") %>%
              select("SampleId", "age"),
            by = "SampleId") %>%
  select(-SampleId) %>%
  group_by(group, age) %>%
  summarise(Influx = sum(RAI, na.rm = TRUE), .groups = "drop")

#herbaceous taxa
herbs_dna <- plant_rai %>%
  filter(group2 %in% c("Novel herbaceous taxa", "Herbaceous taxa", "Nemoral")) %>%
  mutate(group = "Herbs_dna") %>%
  rename_with(~ gsub("RAI_", "", .), starts_with("RAI_")) %>%
  select(group, TRAM_1b_3U_11.12:TRAM_1b_4U_109.110) %>%
  pivot_longer(cols = c("TRAM_1b_3U_11.12":"TRAM_1b_4U_109.110"),
               values_to = 'RAI',
               names_to = 'SampleId') %>%
  left_join(read.csv ("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM21/TRAM21_composite_depth_ages.csv") %>%
              rename("SampleId" = "SampleID",
                     "age" = "median") %>%
              select("SampleId", "age"),
            by = "SampleId") %>%
  group_by(group, age) %>%
  summarise(Influx = sum(RAI, na.rm = TRUE), .groups = "drop")

#steppe taxa
steppe_dna <- plant_rai %>%
  filter(family_name == "Asteraceae" |
           family_name == "Chenopodiaceae" |
           family_name == "Brassicaceae" |
           genus_name == "Helianthemum" |
           genus_name == "Artemisia" |
           genus_name == "Plantago" |
           genus_name == "Rumex" |
           genus_name == "Ephedra") %>%
  rename_with(~ gsub("RAI_", "", .), starts_with("RAI_")) %>%
  mutate(group = "Steppe_dna") %>%
  select(group, TRAM_1b_3U_11.12:TRAM_1b_4U_109.110) %>%
  pivot_longer(cols = c("TRAM_1b_3U_11.12":"TRAM_1b_4U_109.110"),
               values_to = 'RAI',
               names_to = 'SampleId') %>%
  left_join(read.csv ("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM21/TRAM21_composite_depth_ages.csv") %>%
              rename("SampleId" = "SampleID",
                     "age" = "median") %>%
              select("SampleId", "age"),
            by = "SampleId") %>%
  group_by(group, age) %>%
  summarise(Influx = sum(RAI, na.rm = TRUE), .groups = "drop")

#poaceae
poaceae_dna <- plant_rai %>%
  filter(family_name == "Poaceae") %>%
  rename_with(~ gsub("RAI_", "", .), starts_with("RAI_")) %>%
  mutate(group = "Poaceae_dna") %>%
  select(group, TRAM_1b_3U_11.12:TRAM_1b_4U_109.110) %>%
  pivot_longer(cols = c("TRAM_1b_3U_11.12":"TRAM_1b_4U_109.110"),
               values_to = 'RAI',
               names_to = 'SampleId') %>%
  left_join(read.csv ("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM21/TRAM21_composite_depth_ages.csv") %>%
              rename("SampleId" = "SampleID",
                     "age" = "median") %>%
              select("SampleId", "age"),
            by = "SampleId") %>%
  group_by(group, age) %>%
  summarise(Influx = sum(RAI, na.rm = TRUE), .groups = "drop")

#anthropogenic taxa
anthropogenic_dna <- plant_rai %>%
  filter(family_name == "Brassicaceae" |
           species_name == "Urtica dioica" |
           genus_name == "Hordeum" |
           genus_name == "Plantago" |
           genus_name == "Rumex") %>%
  rename_with(~ gsub("RAI_", "", .), starts_with("RAI_")) %>%
  mutate(group = "Anthropogenic_dna") %>%
  select(group, TRAM_1b_3U_11.12:TRAM_1b_4U_109.110) %>%
  pivot_longer(cols = c("TRAM_1b_3U_11.12":"TRAM_1b_4U_109.110"),
               values_to = 'RAI',
               names_to = 'SampleId') %>%
  left_join(read.csv ("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM21/TRAM21_composite_depth_ages.csv") %>%
              rename("SampleId" = "SampleID",
                     "age" = "median") %>%
              select("SampleId", "age"),
            by = "SampleId") %>%
  group_by(group, age) %>%
  summarise(Influx = sum(RAI, na.rm = TRUE), .groups = "drop")

#now i join all herbs DNA into one dataframe
#CAREFUL: the herbaceous component includes poaceae, anthropogenic and steppe, so i will have to exclude them from the herbs_dna dataframe for the figure
herbs_poaceae_steppe_anthropogenic_dna <- herbs_dna %>%
  full_join(poaceae_dna, by = c("group", "age", "Influx")) %>%
  full_join(steppe_dna, by = c("group", "age", "Influx")) %>%
  full_join(anthropogenic_dna, by = c("group", "age", "Influx")) %>%
  pivot_wider(names_from = group, values_from = Influx) %>%
  mutate(Herbs_clean_dna = Herbs_dna - Anthropogenic_dna - Steppe_dna - Poaceae_dna) %>%
  pivot_longer(cols = c(Poaceae_dna, Steppe_dna, Anthropogenic_dna, Herbs_clean_dna),
               names_to = "group",
               values_to = "Influx")

#### 11.1.2. Animal sedaDNA ####
RAI.animals <- read.csv("E:/Saco/IJP/1_data/Tramacastilla/sedaDNA_metabarcoding/RAI/TRAM21.animals.taxa.RAI.nopercentage.csv") %>%
  rename ("SampleId" = "sample_id") %>%
  rename_with(~ gsub("(?<!sp|cf)\\.", " ", .x, perl = TRUE)) %>% #to change only points that are not after "cf" or "sp"
  mutate(across(age:`Nyctalus noctula`, as.numeric))

##### 11.2. Join all data #####
figure_all <- tram20_PAR_groups_discussion %>%
  select(-depth) %>%
  full_join(RAI.animals %>%
              rename("Cattle" = "Bos taurus",
                     "Sheep" = "Ovis aries",
                     "Goat" = "Capra hircus") %>%
              select (age, Cattle, Sheep, Goat), by = "age") %>%
  full_join(CHAR %>%
              select(age, CHAR.tram), by = "age") %>%
  left_join(Mendukilo_binned, by = "age") %>%
  pivot_longer(cols = c("Abies":"d13C.Mendukilo"),
               values_to = 'Influx',
               names_to = 'group') %>%
  filter(age<15000) %>%
  as.data.frame()

XRF_selection <- tram20_XRF_LECO_MS %>%
  right_join(tram20_depths_ages_XRF %>%
               select(SampleId, age), by = "SampleId") %>%
  relocate(age, .after = "composite_depth") %>%
  select(SampleId, age, Si, TOC)

##### 11.3. Plot ######
p.char <- figure_all %>%
  dplyr::filter(group %in% "CHAR.tram") %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2), fill = group)) +
  ylab("")+
  theme_classic() +
  geom_col(width = 0.10) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  guides(fill = guide_legend (title = "")) +
  scale_fill_manual(values = "black", 
                    labels = "CHAR") +
  scale_y_continuous() +
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.text = element_text(size = 11))

p.w.l.ratio <- W.L.ratio %>%
  rename("Type of fuel" = "fuel") %>%
  ggplot(aes(x = age/1000, y = count, fill = `Type of fuel`, color = `Type of fuel`)) +
  geom_area(alpha = 0.6, position = "stack") +
  geom_smooth(method = "loess", se = FALSE, aes(color = `Type of fuel`), linewidth = 1, alpha = 0.8) +
  scale_y_continuous(trans = "sqrt",  # transformation with squared root
                     breaks = c(0, 5, 20, 50)) +
  labs(y = "") +
  theme_classic() +
  scale_x_reverse(limits = c(14.909, min(W.L.ratio$age)/1000),
                  breaks = seq(0.5, 15, by = 0.5),
                  labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("Woody" = "darkgreen", "Grasslands" = "#FFD700")) +
  scale_color_manual(values = c("Woody" = "darkgreen", "Grasslands" = "#FFD700")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.text = element_text(size = 11)) +
  annotate("text", x = 2, y = 120, label = "peak = 59", vjust = -1, color = "red")

p.trees <- figure_all %>%
  dplyr::filter(group %in% c("Abies", "Fagus", "Pinus")) %>%
  mutate(group = fct_relevel(group, "Pinus", "Abies", "Fagus")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "Plant group")) +
  ylab("")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#422CB2", "#BFB2FF", "#E18727")) +
  scale_y_continuous(limits = c(0, 1.25e07),
                     breaks = c(0, 4e+6, 8e+6, 1.2e+7),
                     labels = c("0", "4e+6", "8e+6", "1.2e+7")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", legend.text = element_text(face = "italic", size = 11),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.title = element_text(size = 14))

p.trees_dna <- trees_dna %>%
  dplyr::filter(group %in% c("Abies alba_dna", "Fagus sylvatica_dna", "Pinus_dna")) %>%
  mutate(group = fct_relevel(group, "Pinus_dna", "Abies alba_dna", "Fagus sylvatica_dna")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  ylab("")+
  theme_classic() +
  geom_col_pattern(width = 0.15, pattern = "stripe", pattern_density = 0.15, pattern_spacing = 0.05,
                   pattern_angle = 45, pattern_fill = "black", pattern_colour = NA, aes(pattern_type = "sedaDNA"),
                   show.legend = c(fill = FALSE, pattern = TRUE)) +
  scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#422CB2", "#BFB2FF", "#E18727")) +
  scale_pattern_type_manual(values = "stripe", name = NULL, labels = "sedaDNA") +
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c("0", "0.2", "0.4", "0.6")) +
  theme(axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.position = "none")

p.dec.trees <- figure_all %>%
  dplyr::filter(group %in% c("Betula", "Corylus", "Deciduous Quercus")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "")) +
  ylab("")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#749B58", "#99CC00", "#AE1F63")) +
  scale_y_continuous(limits = c(0, 1.25e07),
                     breaks = c(0, 4e+6, 8e+6, 1.2e+7),
                     labels = c("0", "4e+6", "8e+6", "1.2e+7")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", legend.text = element_text(face = "italic", size = 11),
        plot.margin = margin(0, 0, 0, 0, "cm"))

p.dec.trees_dna <- trees_dna %>%
  dplyr::filter(group %in% c("Betula_dna", "Betulaceae_dna", "Quercus 1_dna")) %>%
  mutate(group = fct_relevel(group, "Betula_dna", "Betulaceae_dna", "Quercus 1_dna")) %>%
  ggplot(aes(x = age/1000, y = Influx, group = group, fill = group)) +
  ylab ("") +
  theme_classic() +
  geom_col_pattern(width = 0.15, pattern = "stripe", pattern_density = 0.15, pattern_spacing = 0.05,
                   pattern_angle = 45, pattern_fill = "black", pattern_colour = NA, aes(pattern_type = "sedaDNA"),
                   show.legend = c(fill = FALSE, pattern = TRUE)) +
  scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                  breaks = seq(0.5, 15, by = 0.5),
                  labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#749B58", "#99CC00", "#AE1F63"), guide = "none") +
  scale_pattern_type_manual(values = "stripe", name = NULL, labels = "sedaDNA") +
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c("0", "0.2", "0.4", "0.6")) +
  theme(axis.title.x = element_blank(),
        legend.position = "right", legend.text = element_text(face = "italic", size = 11),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  guides(pattern_type = guide_legend(title = NULL, override.aes = list(fill = "gray90", pattern = "stripe")))

p.herbs_anthropogenic <- figure_all %>%
  dplyr::filter(group %in% c("Poaceae", "Steppe taxa", "Other herbs", "Anthropogenic taxa")) %>%
  mutate(group = fct_relevel(group, "Poaceae", "Steppe taxa", "Other herbs", "Anthropogenic taxa")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "")) +
  ylab("")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#FFB547", "#FFFF00", "#9ECAE1", "#FF86FF")) +
  scale_y_continuous(limits = c(0, 1.25e07),
                     breaks = c(0, 4e+6, 8e+6, 1.2e+7),
                     labels = c("0", "4e+6", "8e+6", "1.2e+7")) +
  theme(axis.title.x = element_blank(), 
        legend.position = "right",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.text = element_text(size = 11))

p.herbs_dna <- herbs_poaceae_steppe_anthropogenic_dna %>%
  mutate(group = fct_relevel(group, "Poaceae_dna", "Steppe_dna", "Herbs_clean_dna", "Anthropogenic_dna")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "")) +
  ylab("")+
  theme_classic() +
  geom_col_pattern(width = 0.15, pattern = "stripe", pattern_density = 0.15, pattern_spacing = 0.05,
                   pattern_angle = 45, pattern_fill = "black", pattern_colour = NA, aes(pattern_type = "sedaDNA"),
                   show.legend = c(fill = FALSE, pattern = TRUE)) +
  scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#FFB547", "#FFFF00", "#9ECAE1", "#FF86FF")) +
  scale_pattern_type_manual(values = "stripe", name = NULL, labels = "sedaDNA") +
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c("0", "0.2", "0.4", "0.6")) +
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"))

p.xrf.selection <- XRF_selection %>%
  filter(!is.na(TOC)) %>%
  ggplot(aes(x = age / 1000)) +
  geom_line(aes(y = Si/50), color = "#FF7F00", linewidth = 0.5) + # Divided by 50
  geom_line(aes(y = TOC * 100), color = "black", linewidth = 0.5) +  # Multiplied by 100
  labs(x = "") +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000),
                     labels = c("0", "2.5e+4", "5e+4", "7.5e+4", "1e+5"),
                     sec.axis = sec_axis(~./100)) +
  scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                  breaks = seq(0.5, 15, by = 0.5),
                  labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y.left = element_blank(),
        axis.title.y.right = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_color_manual(values = c("#FF7F00", "black"))

p.men <- Mendukilo %>%
  ggplot(aes(x=age/1000, y=d13C.Mendukilo))+
  xlab("Age (ka BP)") +
  ylab("")+
  
  theme_classic() +
  geom_line() + scale_y_reverse() +
  scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                  breaks = seq(0.5, 15, by = 0.5),
                  labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.text = element_text(size = 11))

p.char/p.w.l.ratio/p.trees/p.dec.trees/p.herbs_anthropogenic/p.trees_dna/p.dec.trees_dna/p.herbs_dna/p.xrf.selection/p.men +
  plot_layout(heights = c(1, 1, 2, 2, 2, 1, 1, 1, 1, 2), guides = "collect") +
  theme(plot.margin = margin(0, 0, 0, 0, "pt")) #+ #save pdf A4 vertical
  ggsave("E:/Saco/IJP/5_manuscripts/Cap 2-QSR/figures/All files/discussion_vertical_v2.pdf", height = 11.69, width = 8.27)

###### 12. Discussion figure methodological approach ######
figure_percentages <- tram20_pollen_groups_long %>%
  group_by(SampleId, pollen_group_discussion) %>%
  summarise(total_pollen = sum(pollen_counts, na.rm = TRUE)) %>%
  left_join(tram20_pollen_groups_long %>%
              select(SampleId, depth, age, Lycopodium, `Terrestrial pollen sum`) %>%
              distinct(), by = "SampleId") %>%
  mutate(percentage = (total_pollen/`Terrestrial pollen sum`)*100) %>%
  select(SampleId, depth, age, Lycopodium, pollen_group_discussion, total_pollen, `Terrestrial pollen sum`, percentage) %>%
  rename("group" = "pollen_group_discussion") %>%
  filter(age<15000) %>%
  as.data.frame()

p8.trees_perc <- figure_percentages %>%
  dplyr::filter(group %in% c("Abies", "Fagus", "Pinus")) %>%
  mutate(group = fct_relevel(group, "Pinus", "Abies", "Fagus")) %>%
  ggplot(aes(x=age/1000, y=round(percentage, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "Plant group")) +
  ylab("a) Pollen \n percentages")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#422CB2", "#BFB2FF", "#E18727")) +
  scale_y_continuous(limits = c(0, 75), breaks = c(0, 25, 50, 75)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 14),
        legend.position = "top", legend.text = element_text(face = "italic", size = 11),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.title = element_text(size = 14))

p8.trees_PAR <- figure_all %>%
  dplyr::filter(group %in% c("Abies", "Fagus", "Pinus")) %>%
  mutate(group = fct_relevel(group, "Pinus", "Abies", "Fagus")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "Plant group")) +
  ylab("b) Pollen \n Accumulation \n Rate")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#422CB2", "#BFB2FF", "#E18727")) +
  scale_y_continuous(limits = c(0, 1.25e07),
                     breaks = c(0, 4e+6, 8e+6, 1.2e+7),
                     labels = c("0", "4e+6", "8e+6", "1.2e+7")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 14),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.title = element_text(size = 14))

p8.trees_dna <- trees_dna %>%
  dplyr::filter(group %in% c("Abies alba_dna", "Fagus sylvatica_dna", "Pinus_dna")) %>%
  mutate(group = fct_relevel(group, "Pinus_dna", "Abies alba_dna", "Fagus sylvatica_dna")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  ylab("c) Relative \n Abundance \n Index \n -sedaDNA-") + xlab("Age (ka BP)") +
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#422CB2", "#BFB2FF", "#E18727")) +
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c("0", "0.2", "0.4", "0.6")) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 14),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.position = "none")

p8.dec.trees_perc <- figure_percentages %>%
  dplyr::filter(group %in% c("Betula", "Corylus", "Deciduous Quercus")) %>%
  ggplot(aes(x=age/1000, y=round(percentage, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "")) +
  ylab("")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#749B58", "#99CC00", "#AE1F63")) +
  scale_y_continuous(limits = c(0, 75), breaks = c(0, 25, 50, 75)) +
  theme(axis.title.x = element_blank(),
        legend.position = "top", legend.text = element_text(face = "italic", size = 11),
        plot.margin = margin(0, 0, 0, 0, "cm"))

p8.dec.trees_PAR <- figure_all %>%
  dplyr::filter(group %in% c("Betula", "Corylus", "Deciduous Quercus")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "")) +
  ylab("")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#749B58", "#99CC00", "#AE1F63")) +
  scale_y_continuous(limits = c(0, 1.25e07),
                     breaks = c(0, 4e+6, 8e+6, 1.2e+7),
                     labels = c("0", "4e+6", "8e+6", "1.2e+7")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"))

p8.dec.trees_dna <- trees_dna %>%
  dplyr::filter(group %in% c("Betula_dna", "Betulaceae_dna", "Quercus 1_dna")) %>%
  mutate(group = fct_relevel(group, "Betula_dna", "Betulaceae_dna", "Quercus 1_dna")) %>%
  ggplot(aes(x = age/1000, y = Influx, group = group, fill = group)) +
  ylab ("") + xlab("Age (ka BP)") +
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#749B58", "#99CC00", "#AE1F63"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c("0", "0.2", "0.4", "0.6")) +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"))

p8.herbs_perc <- figure_percentages %>%
  dplyr::filter(group %in% c("Poaceae", "Steppe taxa", "Other herbs", "Anthropogenic taxa")) %>%
  mutate(group = fct_relevel(group, "Poaceae", "Steppe taxa", "Other herbs", "Anthropogenic taxa")) %>%
  ggplot(aes(x=age/1000, y=round(percentage, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "")) +
  ylab("")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#FFB547", "#FFFF00", "#9ECAE1", "#FF86FF")) +
  scale_y_continuous(limits = c(0, 75), breaks = c(0, 25, 50, 75)) +
  theme(axis.title.x = element_blank(), 
        legend.position = "top",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.text = element_text(size = 11))

p8.herbs_PAR <- figure_all %>%
  dplyr::filter(group %in% c("Poaceae", "Steppe taxa", "Other herbs", "Anthropogenic taxa")) %>%
  mutate(group = fct_relevel(group, "Poaceae", "Steppe taxa", "Other herbs", "Anthropogenic taxa")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "")) +
  ylab("")+
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#FFB547", "#FFFF00", "#9ECAE1", "#FF86FF")) +
  scale_y_continuous(limits = c(0, 1.25e07),
                     breaks = c(0, 4e+6, 8e+6, 1.2e+7),
                     labels = c("0", "4e+6", "8e+6", "1.2e+7")) +
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"))

p8.herbs_dna <- herbs_poaceae_steppe_anthropogenic_dna %>%
  mutate(group = fct_relevel(group, "Poaceae_dna", "Steppe_dna", "Herbs_clean_dna", "Anthropogenic_dna")) %>%
  ggplot(aes(x=age/1000, y=round(Influx, 2) , group=group, fill=group))+
  guides (fill = guide_legend(title = "")) +
  ylab("") + xlab("Age (ka BP)") +
  theme_classic() +
  geom_col(width = 0.15) + scale_x_reverse(limits = c(max(figure_all$age)/1000, min(figure_all$age)/1000),
                                           breaks = seq(0.5, 15, by = 0.5),
                                           labels = ifelse(seq(0.5, 15, by = 0.5) %% 1 == 0, as.character(seq(0.5, 15, by = 0.5)), "")) +
  scale_fill_manual(values = c("#FFB547", "#FFFF00", "#9ECAE1", "#FF86FF")) +
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c("0", "0.2", "0.4", "0.6")) +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"))

p8.trees_perc/p8.dec.trees_perc/p8.herbs_perc/p8.trees_PAR/p8.dec.trees_PAR/p8.herbs_PAR/p8.trees_dna/p8.dec.trees_dna/p8.herbs_dna +
  plot_layout(heights = c(3, 4, 4), ncol = 3, nrow = 3) +
  theme(plot.margin = margin(0, 0, 0, 0, "pt")) + #save pdf A4 horizontal
  ggsave("E:/Saco/IJP/5_manuscripts/Cap 2-QSR/figures/All files/discussion_methods_v2.pdf", height = 8.27, width = 11.29)

###### 13. Discussion figure Pyrenean records ######
##### 13.1. Import data #####
#### 13.1.1. Pllan dEstan ####
PDE_pollen <- read_xlsx("E:/Saco/IJP/1_data/Plan d'Estan/PlandEstan_polen.xlsx") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(-c(AP, NAP:PTERIDOPHYTA)) %>%
  mutate(tree = rowSums(select(., Abies:`Quercus perennifolio`), na.rm = TRUE),
         herb = rowSums(select(., Poaceae:Pedicularis), na.rm = TRUE),
         deciduous_trees = rowSums(select(., c(Betula, Corylus, Alnus, Fagus, Tilia, Ulmus, Juglans, Salix, `Quercus caducifolio`,
                                               `Cornus?`)), na.rm = TRUE),
         other_asteraceae = rowSums(select(., c(Cichorioideae, Tubuliflorae, Carduae, Centaurea))),
         steppe = rowSums(select(., c(Myrica, `Ephedra fragilis`, `Ephedra dystachia`, `Cichorioideae`, Tubuliflorae, Carduae, Centaurea,
                                       Artemisia, Chenopodiaceae, Caryophyllaceae, Plantago, Brassicaceae, Rumex, Helianthemum)))) %>%
  pivot_longer(Abies:steppe,
               names_to = "taxon",
               values_to = "percentage")

#### 13.1.2. Marbore ####
marbore_PAR <- read.csv("E:/Saco/IJP/1_data/Marboré/Leunda et al 2019/data/marbore_char_par.csv", sep = ",") %>%
  mutate(tree = rowSums(select(., Abies:Sorbus), na.rm = TRUE),
         herb = rowSums(select(., Poaceae:Typha), na.rm = TRUE)) %>%
  select(cal.BP, Abies:Sorbus, tree, herb) %>%
  pivot_longer(Abies:herb,
               names_to = "taxon",
               values_to = "PAR") %>%
  rename("age" = "cal.BP")

marbore_pollen <- read.csv("E:/Saco/IJP/1_data/Marboré/Leunda et al 2019/data/marbore.csv", sep = ",") %>%
  mutate(pollen_sum = rowSums(select(., Abies:Typha), na.rm = TRUE)) %>%
  rename("age" = "cal.BP") %>%
  select(age, Abies:Typha, pollen_sum) %>%
  mutate(tree = rowSums(select(., Abies:Sorbus), na.rm = TRUE),
         herb = rowSums(select(., Poaceae:Typha), na.rm = TRUE)) %>%
  relocate(pollen_sum, .after = "age") %>%
  pivot_longer(cols = c("Abies":"herb"),
               values_to = "pollen_counts",
               names_to = "taxon") %>%
  mutate_at("pollen_counts", ~replace(., is.na(.), 0)) %>%
  group_by(age, taxon) %>%
  mutate(percentage = (pollen_counts/pollen_sum)*100)

#### 13.1.3. Basa de la Mora ####
bsm_PAR <- read.csv("E:/Saco/IJP/1_data/Marboré/Leunda et al 2019/data/basa_char_par.csv", sep = ",") %>%
  mutate(tree = rowSums(select(., Abies:Pistacia), na.rm = TRUE)) %>%
  select(cal.BP, Abies:Pistacia, tree, herb) %>%
  pivot_longer(Abies:herb,
               names_to = "taxon",
               values_to = "PAR") %>%
  rename("age" = "cal.BP")

bsm_pollen <- read.csv("E:/Saco/IJP/1_data/Marboré/Leunda et al 2019/data/basa.csv", sep = ",") %>%
  mutate(pollen_sum = rowSums(select(., Abies:Pedicularis), na.rm = TRUE)) %>%
  rename("age" = "cal.BP") %>%
  select(age, Abies:Pedicularis, pollen_sum) %>%
  mutate(tree = rowSums(select(., Abies:Pistacia), na.rm = TRUE),
         herb = rowSums(select(., Poaceae:Pedicularis), na.rm = TRUE)) %>%
  relocate(pollen_sum, .after = "age") %>%
  pivot_longer(cols = c("Abies":"herb"),
               values_to = "pollen_counts",
               names_to = "taxon") %>%
  mutate_at("pollen_counts", ~replace(., is.na(.), 0)) %>%
  group_by(age, taxon) %>%
  mutate(percentage = (pollen_counts/pollen_sum)*100)

#### 13.1.4. Bassa Nera ####
bassanera <- read.csv("E:/Saco/IJP/1_data/Bassa Nera/Bassa_Nera_pollen_neotoma.csv", sep = ",") %>%
  select(-c(element, context, units)) %>%
  filter(!group %in% c("TEAM", "ALGA", "FUNG", "AQVP", "VACR", "UNID")) %>% #only keep trees, shrubs and herbs
  pivot_longer(cols = -name, names_to = "SampleId", values_to = "value") %>%
  pivot_wider(names_from = name, values_from = value) %>%
  rename("age" = "Sample age (cal BP)") %>%
  mutate(age = case_when(grepl("--/.+/--", age) ~ gsub("--/(.+)/--", "\\1", age),
                         grepl("--/.+/-", age) ~ gsub("--/(.+)/-", "\\1", age),
                         grepl("-/.+/--", age) ~ gsub("-/(.+)/--", "\\1", age),
                         TRUE ~ age),
         age = gsub("-", "", age),
         age = as.numeric(age)) %>%
  select(SampleId, age, Abies:Ulmus)

bassanera_groups <- bassanera %>%
  filter(SampleId == "group") %>%
  select(-SampleId, -age) %>%
  pivot_longer(everything(), names_to = "taxon", values_to = "group")

bassanera_pollen <-  bassanera %>%
  filter(SampleId != "group") %>%
  select(-SampleId) %>%
  mutate(across(everything(), as.numeric)) %>%
  rowwise() %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(tree = sum(c_across(any_of(bassanera_groups$taxon[bassanera_groups$group == "TR"]))),
         herb = sum(c_across(any_of(bassanera_groups$taxon[bassanera_groups$group == "UPHE"])))) %>%
  ungroup() %>%
  mutate(pollen_sum = rowSums(select(., Abies:Ulmus), na.rm = TRUE)) %>%
  relocate(pollen_sum, .after = "age") %>%
  pivot_longer(cols = c("Abies":"herb"),
               values_to = "pollen_counts",
               names_to = "taxon") %>%
  mutate_at("pollen_counts", ~replace(., is.na(.), 0)) %>%
  group_by(age, taxon) %>%
  mutate(percentage = (pollen_counts/pollen_sum)*100)

#### 13.1.5. Portalet ####
portalet_pollen <- read_csv("E:/Saco/IJP/1_data/Portalet/portalet_final_GGR.csv", col_names = FALSE) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%
  `colnames<-`(.[1,]) %>%
  slice(-1) %>%
  mutate(across(everything(), as.numeric)) %>%
  rename("age" = "Age (cal yr BP)",
         "depth" = "depth (cm)") %>%
  mutate(tree = rowSums(select(., Abies:Pistacia), na.rm = TRUE),
         herb = rowSums(select(., Poaceae:Ranunculaceae), na.rm = TRUE),
         deciduous_trees = rowSums(select(., c(Betula, Corylus, `deciduous Quercus`, Alnus, Salix, Ulmus, Populus, Acer, Fraxinus, Fagus,
                                              Tilia, Juglans, Castanea)), na.rm = TRUE),
         other_asteraceae = rowSums(select(., c(Cichorioideae, Asteroideae, Carduae, Centaurea)), na.rm = TRUE),
         steppe = rowSums(select(., c(Artemisia, Asteroideae, Centaurea, Cichorioideae, Chenopodiaceae,
                                       `Ephedra distachya type`, `Ephedra fragilis type`, Plantago, Helianthemum, Rumex, Caryophyllaceae, Urticaceae)), na.rm = TRUE),
         pollen_sum = rowSums(select(., Abies:Ranunculaceae), na.rm = TRUE)) %>%
  relocate(pollen_sum, .after = "age") %>%
  pivot_longer(cols = c("Abies":"steppe"),
               values_to = "pollen_counts",
               names_to = "taxon") %>%
  mutate_at("pollen_counts", ~replace(., is.na(.), 0)) %>%
  group_by(age, taxon) %>%
  mutate(percentage = (pollen_counts/pollen_sum)*100)

#### 13.1.6. Estanya ####
estanya <- read_xlsx("E:/Saco/IJP/1_data/Estanya/Estanya_polen.xlsx") %>%
  mutate(across(everything(), as.numeric)) %>%
  rename("age" = "Edad",
         "depth" = "Prof. Acumulada",
         "pollen_sum" = "Total Con Pinus, sin acuaticas") %>%
  filter(pollen_sum > 100) %>% #remove samples that have less than 100 pollen grains
  rename_with(~ gsub(",", "", .x)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(-c(Lycopodium, Indeterminado, Isoetes:`Trilete indeterminada`)) %>%
  mutate(across(-c(depth, age, pollen_sum), ~.*100)) #to get percentages, a data is presented by 1

estanya_groups <- estanya %>%
  pivot_longer(Abies:Ruppia, names_to = "taxon", values_to = "pollen_count") %>%
  mutate(group = case_when(taxon %in% c("Betula", "Corylus", "Alnus", "Carpinus", "Salix", "Ulmus",
                                        "Populus", "Acer", "Fagus", "Fraxinus", "Tilia", "Juglans", "Castanea",
                                        "Quercus caducifolio") ~ "deciduous_trees",
                           taxon %in% c("Poaceae", "Artemisia", "Chenopodiaceae", "Caryophyllaceae", "Urticaceae",
                                        "Rumex", "Cichoroideae", "Asteroideae", "Carduae", "Plantago") ~ "steppe",
                           taxon %in% c("Cichorioideae", "Asteroideae", "Carduae", "Centaurea") ~ "other_asteraceae",
                           TRUE ~ "Other"))

estanya_pollen <- estanya %>%
  rowwise() %>%
  mutate(deciduous_trees = sum(c_across(any_of(estanya_groups$taxon[estanya_groups$group == "deciduous_trees"]))),
         steppe = sum(c_across(any_of(estanya_groups$taxon[estanya_groups$group == "steppe"]))),
         other_asteraceae = sum(c_across(any_of(estanya_groups$taxon[estanya_groups$group == "other_asteraceae"])))) %>%
  ungroup() %>%
  relocate(pollen_sum, .after = "age") %>%
  pivot_longer(cols = c("Abies":"other_asteraceae"),
               values_to = "percentage",
               names_to = "taxon")
  
#### 13.1.7. Bosc dels Estanyons ####
estanyons <- read.csv("E:/Saco/IJP/1_data/Bosc dels Estanyons/Estanyons_pollen_neotoma.csv", sep = ",") %>%
  select(-c(element, context, units)) %>%
  filter(!group %in% c("TEAM", "ALGA", "FUNG", "AQVP", "VACR", "UNID", "AQBR", "PTER")) %>% #only keep trees, shrubs and herbs
  pivot_longer(cols = -name, names_to = "SampleId", values_to = "value") %>%
  pivot_wider(names_from = name, values_from = value) %>%
  rename("age" = "author preferred (CAL BP)") %>%
  mutate(age = case_when(grepl("--/.+/--", age) ~ gsub("--/(.+)/--", "\\1", age),
                         grepl("--/.+/-", age) ~ gsub("--/(.+)/-", "\\1", age),
                         grepl("-/.+/--", age) ~ gsub("-/(.+)/--", "\\1", age),
                         TRUE ~ age),
         age = gsub("-", "", age),
         age = as.numeric(age)) %>%
  select(SampleId, age, Abies:Vitis)

estanyons_groups <- estanyons %>%
  filter(SampleId == "group") %>%
  select(-SampleId, -age) %>%
  pivot_longer(everything(), names_to = "taxon", values_to = "group") %>%
  mutate(group = case_when(taxon %in% c("Abies", "Acer", "Alnus", "Betula", "Carpinus", "Castanea", "Corylus",
                                        "Fagus", "Fraxinus", "Juglans", "Olea", "Ostrya", "Pinus", "Pistacia",
                                        "Populus", "Quercus", "Quercus ilex", "Salix", "Sorbus", "Tilia", "Ulmus") ~ "TR",
                           taxon %in% c("Buxus", "Calluna", "Cistus", "Daphne", "Ephedra", "Erica", "Ericaceae",
                                        "Genista/Cytisus", "Hedera", "Ilex aquifolium", "Juniperus", "Lonicera",
                                        "Phillyrea", "Prunus", "Rhamnus", "Rhododendron", "Sambucus", "Vaccinium",
                                        "Viburnum", "Viscum", "Vitis") ~ "SH",
                           group == "UPHE" ~ "UPHE"))

estanyons_pollen <- estanyons %>%
  filter(SampleId != "group") %>%
  select(-SampleId) %>%
  mutate(across(everything(), as.numeric)) %>%
  rowwise() %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(tree = sum(c_across(any_of(estanyons_groups$taxon[estanyons_groups$group == "TR"]))),
         herb = sum(c_across(any_of(estanyons_groups$taxon[estanyons_groups$group == "UPHE"])))) %>%
  ungroup() %>%
  mutate(pollen_sum = rowSums(select(., Abies:Ulmus), na.rm = TRUE)) %>%
  relocate(pollen_sum, .after = "age") %>%
  pivot_longer(cols = c("Abies":"herb"),
               values_to = "pollen_counts",
               names_to = "taxon") %>%
  mutate_at("pollen_counts", ~replace(., is.na(.), 0)) %>%
  group_by(age, taxon) %>%
  mutate(percentage = (pollen_counts/pollen_sum)*100)

##### 13.2. Plot Holocene figure #####
periods_lines_holocene <- c(4.2, 8.2, 11.7) #to add Holocene periods
p11.tramacastilla.v <- tram20_pollen %>%
  select(-c(`No of taxa`, Lycopodium, `Pteropsida trilete undiff`:`Charcoal <150`)) %>%
  pivot_longer(cols = c("Abies":"Hygrophytes"),
               values_to = 'pollen_counts',
               names_to = 'taxon') %>%
  mutate_at("pollen_counts", ~replace(., is.na(.), 0)) %>%
  group_by(SampleId, taxon) %>%
  mutate(percentage = (pollen_counts/`Terrestrial pollen sum`)*100) %>%
  filter(taxon %in% c("Betula", "Corylus", "Deciduous Quercus", "Abies", "Pinus", "Herbs")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Abies", "Pinus", "Betula", "Corylus", "Deciduous Quercus", "Herbs")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(11.7, 0),
                  breaks = seq(0, 12, by = 0.5),
                  labels = ifelse(seq(0, 12, by = 0.5) %% 1 == 0, as.character(seq(0, 12, by = 0.5)), "")) +
  scale_fill_manual(values = c("#BFB2FF", "#422CB2", "#749B58", "#99CC00", "#AE1F63", "#9ECAE1")) +
  theme_classic() +
  labs(x = "Age (ka BP)", title = "Tramacastilla \n (1682 m a.s.l.)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_holocene, color = "black", linetype = "dashed", linewidth = 0.3)

p11.portalet.v <- portalet_pollen %>%
  filter(taxon %in% c("Betula", "Corylus", "deciduous Quercus", "Abies", "Pinus", "herb")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Abies", "Pinus", "Betula", "Corylus", "deciduous Quercus", "herb")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(11.7, 0),
                  breaks = seq(0, 12, by = 0.5),
                  labels = ifelse(seq(0, 12, by = 0.5) %% 1 == 0, as.character(seq(0, 12, by = 0.5)), "")) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#BFB2FF", "#422CB2", "#749B58", "#99CC00", "#AE1F63", "#9ECAE1")) +
  theme_classic() +
  labs(title = "El Portalet peatbog \n (1850 m a.s.l.)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_holocene, color = "black", linetype = "dashed", linewidth = 0.3)

p11.marbore.v <- marbore_pollen %>%
  filter(taxon %in% c("Betula", "Corylus", "Dec_Querc", "Abies", "Pinus", "herb")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Abies", "Pinus", "Betula", "Corylus", "Dec_Querc", "herb")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(11.7, 0),
                  breaks = seq(0, 12, by = 0.5),
                  labels = ifelse(seq(0, 12, by = 0.5) %% 1 == 0, as.character(seq(0, 12, by = 0.5)), "")) +
  scale_fill_manual(values = c("#BFB2FF", "#422CB2", "#749B58", "#99CC00", "#AE1F63", "#9ECAE1")) +
  theme_classic() +
  labs(title = "Marboré \n (2612 m a.s.l.)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_holocene, color = "black", linetype = "dashed", linewidth = 0.3)

p11.estanyons.v <- estanyons_pollen %>%
  filter(taxon %in% c("Betula", "Corylus", "Quercus", "Abies", "Pinus", "herb")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Abies", "Pinus", "Betula", "Corylus", "Quercus", "herb")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(11.7, 0),
                  breaks = seq(0, 12, by = 0.5),
                  labels = ifelse(seq(0, 12, by = 0.5) %% 1 == 0, as.character(seq(0, 12, by = 0.5)), "")) +
  scale_fill_manual(values = c("#BFB2FF", "#422CB2", "#749B58", "#99CC00", "#AE1F63", "#9ECAE1")) +
  theme_classic() +
  labs(title = "Bosc dels Estanyons \n (2180 m a.s.l.)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_holocene, color = "black", linetype = "dashed", linewidth = 0.3)

p11.bsm.v <- bsm_pollen %>%
  filter(taxon %in% c("Betula", "Corylus", "Dec_Querc", "Abies", "Pinus", "herb")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Abies", "Pinus", "Betula", "Corylus", "Dec_Querc", "herb")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(11.7, 0),
                  breaks = seq(0, 12, by = 0.5),
                  labels = ifelse(seq(0, 12, by = 0.5) %% 1 == 0, as.character(seq(0, 12, by = 0.5)), "")) +
  scale_fill_manual(values = c("#BFB2FF", "#422CB2", "#749B58", "#99CC00", "#AE1F63", "#9ECAE1")) +
  theme_classic() +
  labs(title = "Basa de la Mora \n (1914 m a.s.l.)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_holocene, color = "black", linetype = "dashed", linewidth = 0.3)

p11.bassanera.v <- bassanera_pollen %>%
  filter(taxon %in% c("Betula", "Corylus", "Quercus (deciduous)", "Abies", "Pinus", "herb")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Abies", "Pinus", "Betula", "Corylus", "Quercus (deciduous)", "herb")) %>%
  mutate(taxon = fct_recode(taxon,
                            "Deciduous~italic(Quercus)" = "Quercus (deciduous)",
                            "italic(Abies)" = "Abies",
                            "italic(Pinus)" = "Pinus",
                            "italic(Betula)" = "Betula",
                            "italic(Corylus)" = "Corylus",
                            "Herbaceous~taxa" = "herb")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(11.7, 0),
                  breaks = seq(0, 12, by = 0.5),
                  labels = ifelse(seq(0, 12, by = 0.5) %% 1 == 0, as.character(seq(0, 12, by = 0.5)), "")) +
  scale_fill_manual(values = c("#BFB2FF", "#422CB2", "#749B58", "#99CC00", "#AE1F63", "#9ECAE1"),
                    labels = c(expression(italic(Abies)), expression(italic(Pinus)), expression(italic(Betula)),
                               expression(italic(Corylus)), expression("Deciduous"~italic(Quercus)), expression("Herbaceous taxa"))) +
  guides (fill = guide_legend(title = "Plant group")) +
  theme_classic() +
  labs(y = "Pollen percentages", title = "Bassa Nera \n (1891 m a.s.l.)") +
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_holocene, color = "black", linetype = "dashed", linewidth = 0.3)

p11.marbore.v/p11.bsm.v/p11.bassanera.v/p11.portalet.v/p11.tramacastilla.v +
  plot_layout(ncol = 1, nrow = 5, guides = "collect") +
  theme(plot.margin = margin(0, 0, 0, 0, "pt")) #+ #save pdf A4 horizontal
  ggsave("E:/Saco/IJP/5_manuscripts/Cap 2-QSR/figures/All files/discussion_records_holocene.pdf", height = 8.27, width = 11.69)

##### 13.3 Plot Lateglacial records #####
periods_lines_lateglacial <- c(11.7, 12.9, 14, 14.9) # to add Lateglacial periods
p12.portalet <- portalet_pollen %>%
  filter(taxon %in% c("deciduous_trees", "Juniperus", "steppe", "Artemisia", "other_asteraceae")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "deciduous_trees", "Juniperus", "steppe", "Artemisia", "other_asteraceae")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(15, 11),
                  breaks = seq(11, 15, by = 0.5),
                  labels = function(x) ifelse(x == floor(x), as.integer(x), x)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#749B58", "#993D00", "#FFFF00", "#7EC3E5", "#FF99BF")) +
  theme_classic() +
  labs(title = "El Portalet peatbog \n (1850 m a.s.l.)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_lateglacial, color = "black", linetype = "dashed", linewidth = 0.3)
  
p12.pde <- PDE_pollen %>%
  filter(taxon %in% c("deciduous_trees", "Juniperus", "steppe", "Artemisia", "other_asteraceae")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "deciduous_trees", "Juniperus", "steppe", "Artemisia", "other_asteraceae")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(15, 11),
                  breaks = seq(11, 15, by = 0.5),
                  labels = function(x) ifelse(x == floor(x), as.integer(x), x)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#749B58", "#993D00", "#FFFF00", "#7EC3E5", "#FF99BF")) +
  theme_classic() +
  labs(title = "Pllan d'Están \n (1840 m a.s.l.)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_lateglacial, color = "black", linetype = "dashed", linewidth = 0.3)

p12.tram <- tram20_pollen %>%
  select(-c(`No of taxa`, Lycopodium, `Pteropsida trilete undiff`:`Charcoal <150`)) %>%
  rowwise() %>%
  mutate(steppe = sum(c_across(c(`Ephedra distachya`, `Ephedra fragilis`, `Asteraceae Artemisia`, `Asteraceae Centaurea`, `Asteraceae Cichorioideae`,
                                `Asteraceae Cichorioideae`, `Asteraceae Carduae`, Brassicaceae, Chenopodiaceae, `Cistaceae t Helianthemum`)), na.rm = TRUE),
         other_asteraceae = sum(c_across(c(`Asteraceae Centaurea`, `Asteraceae Cichorioideae`,
                                         `Asteraceae Cichorioideae`, `Asteraceae Carduae`)), na.rm = TRUE)) %>%
  ungroup() %>%
  rename("deciduous_trees" = "Deciduous trees",
         "Artemisia" = "Asteraceae Artemisia") %>%
  pivot_longer(cols = c("Abies":"other_asteraceae"),
               values_to = 'pollen_counts',
               names_to = 'taxon') %>%
  mutate_at("pollen_counts", ~replace(., is.na(.), 0)) %>%
  group_by(SampleId, taxon) %>%
  mutate(percentage = (pollen_counts/`Terrestrial pollen sum`)*100) %>%
  filter(taxon %in% c("deciduous_trees", "Juniperus", "steppe", "Artemisia", "other_asteraceae")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "deciduous_trees", "Juniperus", "steppe", "other_asteraceae", "Artemisia")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(15, 11),
                  breaks = seq(11, 15, by = 0.5),
                  labels = function(x) ifelse(x == floor(x), as.integer(x), x)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#749B58", "#993D00", "#FFFF00", "#7EC3E5", "#FF99BF"),
                    labels = c(expression("Deciduous trees"), expression(italic(Juniperus)),
                               expression("Steppe taxa"), expression("Other Asteraceae"), expression(italic(Artemisia)))) +
  guides (fill = guide_legend(title = "Plant group")) +
  theme_classic() +
  labs(title = "Tramacastilla lake \n (1682 m a.s.l.)", y = "Pollen percentages") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_lateglacial, color = "black", linetype = "dashed", linewidth = 0.3)

p12.estanya <- estanya_pollen %>%
  filter(taxon %in% c("deciduous_trees", "Juniperus", "steppe", "Artemisia", "other_asteraceae")) %>%
  mutate(taxon = factor(taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "deciduous_trees", "Juniperus", "steppe", "Artemisia", "other_asteraceae")) %>%
  ggplot(aes(x = age/1000, y = percentage, group = taxon, fill = taxon)) +
  geom_area()+
  scale_x_reverse(limits = c(15, 11),
                  breaks = seq(11, 15, by = 0.5),
                  labels = function(x) ifelse(x == floor(x), as.integer(x), x)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#749B58", "#993D00", "#FFFF00", "#7EC3E5", "#FF99BF")) +
  theme_classic() +
  labs(title = "Estanya \n (670 m a.s.l.)", x = "Age (ka BP)") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title.position = "panel") +
  geom_vline(xintercept = periods_lines_lateglacial, color = "black", linetype = "dashed", linewidth = 0.3)

p12.portalet/p12.pde/p12.tram/p12.estanya +
  plot_layout(ncol = 1, nrow = 4, guides = "collect") +
  theme(plot.margin = margin(0, 0, 0, 0, "pt")) + #save pdf A4 horizontal
  ggsave("E:/Saco/IJP/5_manuscripts/Cap 2-QSR/figures/All files/discussion_records_lateglacial.pdf", height = 8.27, width = 11.69)


