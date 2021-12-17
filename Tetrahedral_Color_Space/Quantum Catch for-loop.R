#####Quantum Catch for-loop 

setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/Averaged.Spec.readings")

source("RtcsFunctions.R")	
load("bluetitss.dat")
source("visualization.R")

library(pavo)
library(phytools)
library(tidyverse)
library(pavo)
library(phytools)
library(furrr)
library(future)

Cardinals <- c("Amaur_con_aeq_F","Amaur_con_aeq_M","Amaur_con_rel_F","Amaur_con_rel_M","Amaur_moe_F","Amaur_moe_M","Card_card_F","Card_card_M","Card_pho_F","Card_pho_M","Card_sin_F","Card_sin_M","Cary_cana_fro_F","Cary_cana_fro_M","Cary_pol_F","Cary_pol_M","Chl_car_F","Chl_carr_M","Chl_oli_F","Chl_oli_M","Chl_stolz_F","Chl_stolz_M","Cyan_cyan_F","Cyan_cyan_M","Cyan_glau_F","Cyan_glau_M","Cyan_pare_F","Cyan_pare_M","Cyan_roth_F","Cyan_roth_M","Gran_pel_pel_F","Gran_pel_pel_M","Gran_sal_F","Gran_sal_M","Gran_ven_F","Gran_ven_M","Hab_atrim_F","Hab_atrim_M","Hab_cris_F","Hab_cris_M","Hab_fusc_F","Hab_fusc_M","Hab_gut_F","Hab_gut_M","Hab_rubica_F","Hab_rubica_M","Pass_amo_F","Pass_amo_M","Pass_bris_F","Pass_bris_M","Pass_caer_f","Pass_caer_m","Pass_ciris_F","Pass_ciris_M","Pass_cya_F","Pass_cya_M","Pass_lec_F","Pass_lec_M","Pass_ros_F","Pass_ros_M","Pass_veri_F","Pass_veri_M","Peri_ery_F","Peri_ery_M","Pheu_auro_F","Pheu_auro_M","Pheu_chry_aura_F","Pheu_chry_aura_M","Pheu_chry_F","Pheu_chry_M","Pheu_gas_F","Pheu_gas_M","Pheu_ludo_F","Pheu_ludo_M","Pheu_mel_F","Pheu_mel_M","Pheu_tib_F","Pheu_tib_M","Pir_bid_F","Pir_bid_M","Pir_ery_F","Pir_ery_M","Pir_flava_F","Pir_flava_M","Pir_hep_F","Pir_hep_M","Pir_leu_F","Pir_leu_M","Pir_ludo_F","Pir_ludo_M","Pir_lutea_F","Pir_lutea_M","Pir_oli_F","Pir_oli_M","Pir_rose_F","Pir_rose_M","Pir_rubra_F","Pir_rubra_M","Pir_ruceps_F","Pir_ruceps_M","Rhodo_cel_F","Rhodo_cel_M","Spiz_amer_F","Spiz_amer_M")


# SETTING CORES -- for parrelelizing code ---------------------------------
# sets the number of cores to one fewer than the max available
num_cores <- availableCores() -1
# actually sets the planning algoreithm for how the code will be mapped across the cores.
plan(multiprocess, workers = num_cores) 


# READING -- all CSV files in working directory ---------------------------
files <- dir(pattern = "*.csv")

reads <- files%>%
  purrr::map(read_csv)

overlap_patches <- c("Crown","Auriculars","Nape","Mantle","Back","Rump","Dorsal.tail","Throat","Breast","Side","Belly","Undertail.coverts","Scapulars","Medium.Coverts","Greater.Coverts","Primaries","Secondaires")


vis_model <-reads%>%
  future_map(~ gather(., code, values, 2:ncol(.))%>%
               #mutate(values = abs(values))%>%
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
                                               TRUE ~ as.character(speices_code)))%>%
               filter(patch %in% overlap_patches)%>%
               #dplyr::select(c(wl, overlap_patches))%>%
               separate(speices_code, c("species", "sex"), sep = -1)%>%
               spread(patch, values)%>%
               group_by(species, sex)%>%
               nest())%>%
  reduce(bind_rows)%>%
  mutate(data = future_map(data, ~ as.data.frame(.x)%>%
                             vismodel(., visual = "avg.uv", achromatic = "bt.dc")))%>%
  unnest(data)

write_csv(vis_model, "quatnum.dat")











###### Pavo dichromatism attempt

library(pavo)
#Read in the reflectance spectra
refs1 <- read.csv("Amaur_moe_F_output.csv")
refs1 <- read.csv("Card_card_M_output.csv")
refs <- refs1[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]

vismod1 <- vismodel(refs1,   visual = "avg.uv", achromatic = "bt.dc")
vismod2 <- vismodel(refs2,   visual = "avg.uv", achromatic = "bt.dc")
vismod1

#allrefs <- cbind(vismod1,vismod2)

colddist1 <- coldist(vismod1,
        noise = "neural", achromatic = TRUE, n = c(1, 2, 2, 4),
        weber = 0.1, weber.achro = 0.1)

card <- coldist(data, subset = c('PirLeuM','PirLeuF'))

##### for-loop

filenames <- list.files(path=getwd(),pattern="*.csv") 
numfiles <- length(filenames)  
datalist = list()

for (i in c(1:numfiles)){
  tryCatch({
    print(filenames[i])
    refs <- read.csv(filenames[i], header=TRUE)
    refs <- refs[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]
    vismod1 <- vismodel(refs,   visual = "bluetit", achromatic = "bt.dc", trans = "bluetit")
    datalist[[i]] <- vismod1},
    error = function(err) {
      
      print(paste("Didn't work:  ",filenames[i]))
      
    },
    warning = function(warn) {}, finally = {}
  )
}
data <- do.call(rbind, datalist)
write.csv(data, file = "quatnum.csv")
