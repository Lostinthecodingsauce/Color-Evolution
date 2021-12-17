## Analyzing patch variability using vis_model
## written by Zach Quinlan for Ben Scott
## 10/22/2019
setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/Averaged.Spec.readings")
# READING -- libraries ------------------------------------------------------
library(tidyverse)
library(pavo)
library(phytools)
library(furrr)
library(future)
library(dplyr)
library(rgl)
source("RtcsFunctions.R")	
source("More.Functions.R")
load("bluetitss.dat")


# SETTING CORES -- for parrelelizing code ---------------------------------
# sets the number of cores to one fewer than the max available
num_cores <- availableCores() -1
# actually sets the planning algoreithm for how the code will be mapped across the cores.
plan(multiprocess, workers = num_cores) 


# READING -- all CSV files in working directory ---------------------------
files <- dir(pattern = "*.csv")

reads <- files%>%
  purrr::map(read_csv)

mtree <- read.nexus("Cardinalidae.nexus")

overlap_patches <- c("Crown","Auriculars","Nape","Mantle",
                     "Back","Rump","Dorsal.tail","Throat",
                     "Breast","Side","Belly","Undertail.coverts",
                     "Scapulars","Medium.Coverts","Greater.Coverts"
                     ,"Secondaires")

bad_birds <- c("AmaurConAeq",
               "AmaurMoe",
               "HabCris",
               "PheuChry",
               "PirHep",
               "PirLutea")

bird_names <- read_csv("./second/name_fix.csv")%>%
  dplyr::select(-c(1,4))

# Models ------------------------------------------------------------------
vis_model <-reads%>%
  future_map(~ gather(., code, values, 2:ncol(.))%>%
               mutate(values = abs(values))%>%
               separate(code, c("patch", "speices_code"), sep = "_")%>%
               mutate(patch = case_when(patch == "Ventrail.tail" ~ "Ventral.tail",
                                        patch == "Undertail" ~ "Ventral.tail",
                                        patch == "PirFaries" ~ "Pirmaries",
                                        patch == "Wingbars" ~ "Wingbar",
                                        patch == "Beak.2" ~ "Back.2",
                                        patch == "Fantle" ~ "Mantle",
                                        patch == "FediuF.Coverts" ~ "Medium.Coverts",
                                        patch == "Mante.2" ~ "Mantle.2",
                                        patch == "RuFp" ~ "Rump",
                                        TRUE ~ as.character(patch)),
                      speices_code = case_when(speices_code == "GranSalF" ~ "GranSelF",
                                               speices_code == "GranSalM" ~ "GranSelM",
                                               speices_code == "CardsinM" ~ "CardSinM",
                                               TRUE ~ as.character(speices_code)))%>%
               filter(patch %in% overlap_patches)%>%
               separate(speices_code, c("species", "sex"), sep = -1)%>%
               filter(!species %in% bad_birds)%>%
               group_by(species, sex, patch)%>%
               nest())%>%
  reduce(bind_rows)%>%
  mutate(data = purrr::map(data, ~ as.data.frame(.x)%>%
                                 vismodel(., visual = "avg.uv", achromatic = "bt.dc")))%>%
  unnest(data)

ace_model <- vis_model%>%
  ungroup()%>%
  full_join(., bird_names, by = "species")%>%
  dplyr::select(-c(species, lum))%>%
  rename(species = species_2)%>%
  dplyr::select(species, sex, patch, everything())%>%
  gather(cone, value, 4:ncol(.))%>%
  group_by(sex, patch, cone)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~
                             (as.character(.x$value)%>%
                                     as.numeric()%>%
                                     ace(., phy = mtree, type="continuous", method="REML", model="BM"))$ace%>%
                             as.data.frame()))%>%
  dplyr::select(-data)%>%
  unnest(model)


#write_csv(vis_model, "quatnum.dat")
write_csv(ace_model, "./second/ace_model.dat")
