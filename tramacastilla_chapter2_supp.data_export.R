#Export Supp Data
#Needed version R 4.2.3
load("E:/Saco/IJP/2_code/tramacastilla_chapter2/tramacastilla_multiproxy.RData")
library(tidyverse)
library(openxlsx)
supp_data1 <- pollen_group %>%
  select(Taxon, pollen_group_general, pollen_group_discussion)
supp_data2 <- tram20_pollen %>%
  select(-c(`Evergreen trees`:Trees, Shrubs, Herbs, Undetermined, Hygrophytes, Hygrophytes, Spores, `Charcoal <150`)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  rename( "Depth (cm)" = "depth") %>%
  rename("Age (cal yr BP)" = "age") %>%
  rename("Number of taxa" = "No of taxa")
supp_data3 <- tram20_PAR %>%
  left_join(tram20_depth_ages %>%
              rename ("age" = "median") %>%
              select(SampleId, age),
            by = "age") %>%
  relocate(SampleId, .before = depth) %>%
  rename( "Depth (cm)" = "depth") %>%
  rename("Age (cal yr BP)" = "age")
supp_data4 <- CharTotS.110.age %>%
  select(-SampleId, -depth) %>%
  full_join(CHAR %>% select(age, CHAR.tram), by = "age") %>%
  full_join(tram20_depth_ages_sed_rate %>% select(SampleId, age, depth), by = "age") %>%
  arrange(age) %>%
  relocate(SampleId, .before = "CharCount") %>%
  relocate(age, .after = "SampleId") %>%
  relocate(depth, .after = "age") %>%
  rename("Age (cal yr BP)" = "age") %>%
  rename( "Depth (cm)" = "depth") %>%
  rename("Number of charcoal particles" = "CharCount") %>%
  rename("Sum of charcoal area" = "sumCharArea") %>%
  rename("mean W/L ratio" = "meanCharWL") %>%
  rename("median W/L ratio" = "medianCharWL") %>%
  rename("CHAR (mm2·cm-2·yr-1)" = "CHAR.tram")
supp_data5 <- tram20_diatoms %>%
  rename( "Depth (cm)" = "depth") %>%
  rename("Age (cal yr BP)" = "age") %>%
  rename("Number of taxa" = "richness")
supp_data6 <- trees_dna %>%
  full_join(herbs_poaceae_steppic_anthropogenic_dna %>%
              select(-Herbs_dna),by = c("age", "Influx", "group")) %>%
  arrange(age) %>%
  pivot_wider(names_from = group, values_from = Influx) %>%
  rename_with(~ gsub("_dna", "", .x)) %>%
  rename("Corylus" = "Betulaceae",
         "Deciduous Quercus" = "Quercus 1",
         "Herbaceous taxa" = "Herbs_clean")
supp_data7 <- read.csv("E:/Saco/IJP/1_data/Tramacastilla/sedaDNA_metabarcoding/RAI/TRAM21.animals.groups.RAI.nopercentage.csv") %>%
  rename ("SampleId" = "sample_id") %>%
  left_join(read.csv ("E:/Saco/IJP/1_data/Tramacastilla/depth/TRAM21/TRAM21_composite_depth_ages.csv") %>%
              rename("SampleId" = "SampleID",
                     "age" = "median") %>%
              select("SampleId", "age"),
            by = "SampleId") %>%
  rename_with(~ gsub("d.", "d ", .x)) %>%
  rename_with(~ gsub("r.", "r ", .x)) %>%
  select(-SampleId) %>%
  relocate(age, .before = 1)
supp_data8 <- tram20_depth_ages_sed_rate %>%
  select(-c(depth_diff, age_diff)) %>%
  rename("median" = "age") %>%
  rename("sedimentation rate (yrs/cm)" = "sed_rate")
supp_data9 <- read.csv("E:/Saco/Tramacastilla/Maps/Coordenates_references.csv") %>%
  filter(!name %in% c("Lago della Costa", "Lago Sangiatto")) %>%
  select(-ID) %>%
  rename("Altitude (m a.s.l.)" = "Altitude")

suppdata <- createWorkbook()

addWorksheet(suppdata, "Supp. Data 1-Pollen types class")
writeData(suppdata, "Supp. Data 1-Pollen types class", supp_data1)

addWorksheet(suppdata, "Supp. Data 2-Pollen counts")
writeData(suppdata, "Supp. Data 2-Pollen counts", supp_data2)

addWorksheet(suppdata, "Supp. Data 3-PAR values")
writeData(suppdata, "Supp. Data 3-PAR values", supp_data3)

addWorksheet(suppdata, "Supp. Data 4-Charcoal metrics")
writeData(suppdata, "Supp. Data 4-Charcoal metrics", supp_data4)

addWorksheet(suppdata, "Supp. Data 5-Diatom counts")
writeData(suppdata, "Supp. Data 5-Diatom counts", supp_data5)

addWorksheet(suppdata, "Supp. Data 6-Plant sedaDNA")
writeData(suppdata, "Supp. Data 6-Plant sedaDNA", supp_data6)

addWorksheet(suppdata, "Supp. Data 7-Animal sedaDNA")
writeData(suppdata, "Supp. Data 7-Animal sedaDNA", supp_data7)

addWorksheet(suppdata, "Supp. Data 8-Ages and depths")
writeData(suppdata, "Supp. Data 8-Ages and depths", supp_data8)

addWorksheet(suppdata, "Supp. Data 9-Coordinates")
writeData(suppdata, "Supp. Data 9-Coordinates", supp_data9)

saveWorkbook(suppdata, "E:/Saco/IJP/5_manuscripts/Cap 2-QSR/supplementary material/Supplementary Data.xlsx",
             overwrite = TRUE)

