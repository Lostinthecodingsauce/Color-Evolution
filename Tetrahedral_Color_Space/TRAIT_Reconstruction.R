######   ANCESTRAL STATE RECONSTRUCTION  ##########################

setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis")
source("Code/PLumage.ancestral.recon.R")	
source("Code/visualization.R")
load("bluetitss.dat")

# Read in male Quantum Cone catch dataframe, read in female Quantum Cone catch dataframe 
all_male_qcatch <- read.csv("m_qcatch.csv", header = TRUE)
all_female_qcatch <- read.csv("f_qcatch.csv", header = TRUE)
#Read in Tree
tree <- read.nexus("Cardinalidae.nexus")




###### need to subset each plumage region for analysis. Yes, each fucking region
patches <- c("Crown","Auriculars","Nape","Mantle","Back","Rump","Dorsal.tail","Throat","Breast","Side","Belly","Undertail.coverts","Scapulars","Medium.Coverts","Greater.Coverts","Primaries","Secondaires")

# Test of single region
m_crown <- subset(all_male_qcatch, Patch == "Crown")
rownames(m_crown) <- m_crown[,1,]
m_crown <- m_crown[tree$tip.label, ]  

f_crown <- subset(all_female_qcatch, Patch == "Crown")
rownames(f_crown) <- f_crown[,1,]
muir_region_mf <- as.data.frame(rbind(m_crown,f_crown)[,1:7])


price_eaton_asr(tree, m_crown, f_crown)


#male region
muir_u_asr_m <- ace(phy=tree,x=as.numeric(as.character(m_crown$u)),type="continuous",method="REML",model="BM")
muir_s_asr_m <- ace(phy=tree,x=as.numeric(as.character(m_crown$s)),type="continuous",method="REML",model="BM")
muir_m_asr_m <- ace(phy=tree,x=as.numeric(as.character(m_crown$m)),type="continuous",method="REML",model="BM")
muir_l_asr_m <- ace(phy=tree,x=as.numeric(as.character(m_crown$l)),type="continuous",method="REML",model="BM")
#female region 
muir_u_asr_f <- ace(phy=tree,x=as.numeric(as.character(f_crown$u)),type="continuous",method="REML",model="BM")
muir_s_asr_f <- ace(phy=tree,x=as.numeric(as.character(f_crown$s)),type="continuous",method="REML",model="BM")
muir_m_asr_f <- ace(phy=tree,x=as.numeric(as.character(f_crown$m)),type="continuous",method="REML",model="BM")
muir_l_asr_f <- ace(phy=tree,x=as.numeric(as.character(f_crown$l)),type="continuous",method="REML",model="BM")


#Combine all ancestral states into a dataframe so we can calculate color distances
muir_node_qc_m <- data.frame(muir_u_asr_m$ace,muir_s_asr_m$ace,muir_m_asr_m$ace,muir_l_asr_m$ace)
rownames(muir_node_qc_m) <- paste(names(muir_u_asr_m$ace),"_m",sep="")
rownames(muir_node_qc_m) <- rownames(m_crown)=rownames(muir_node_qc_m)

colnames(muir_node_qc_m) <- c("u","s","m","l")
muir_node_qc_f <- data.frame(muir_u_asr_f$ace,muir_s_asr_f$ace,muir_m_asr_f$ace,muir_l_asr_f$ace)
rownames(muir_node_qc_f) <- paste(names(muir_u_asr_m$ace),"_f",sep="")
colnames(muir_node_qc_f) <- c("u","s","m","l")
muir_node_qc_mf <- as.data.frame(rbind(muir_node_qc_m,muir_node_qc_f))

names_muir <- rownames(muir_region_mf)
muir_region_mf <- apply(muir_region_mf,2,as.numeric)
rownames(muir_region_mf) <- names_muir



muir_node_coldist <- coldist(muir_node_qc_mf,qcatch="Qi",vis="tetra",noise="neural",subset=NULL,achro=T,n1=1,n2=2,n3=2,n4=4,v=0.05)

muir_node_coldist <- coldist(muir_node_qc_mf,qcatch="Qi",noise="neural",subset=NULL,
                             achro=T,n = c(1,2,3,4,5,6,7,8,9,10,11),weber = 0.1, weber.achro = 0.1)


muir_node_coldist



patches <- c("Crown","Auriculars","Nape","Mantle","Back","Rump","Dorsal.tail","Throat","Breast","Side","Belly","Undertail.coverts","Scapulars","Medium.Coverts","Greater.Coverts","Primaries","Secondaires")


filter(patch %in% overlap_patches)%>%

# Test of single region
m_crown <- subset(all_male_qcatch, ï..Patch == "Crown")
f_crown <- subset(all_female_qcatch, ï..Patch == "Crown")

datalist = list()

for (i in c(1:all_male_qcatch)){
  tryCatch({
    print(all_male_qcatch[i]))
  
  datalist[[i]] <- },
  error = function(err) {
    
    print(paste("Didn't work:  ",filenames[i]))
    
  },
  warning = function(warn) {}, finally = {}
  )
}


### Comebine results and write them as csv
data <- do.call(rbind, datalist)


filter(i..patch %in% overlap_patches)%>%
  
    trivalues <- tristimulus(refs, ss)
    datalist[[i]] <- trivalues}
 





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
                              as.data.frame()%>%
                              rownames_to_column(var = "node_number")%>%
                              rename("value" = 2)))%>%
  dplyr::select(-data)%>%
  unnest(model)%>%
  spread(cone, value)

col_dist <-ace_model%>%
  group_by(sex, patch)%>%
  nest()%>%
  mutate(data = purrr::map(data, ~ 
                             column_to_rownames(.x, var =  "node_number")%>%
                             coldist(qcatch="Qi", noise="neural", subset=NULL, achro=T, n = c(1,2,2,4), weber = 0.05)))

write_csv(vis_model, "quatnum.dat")
write_csv(ace_model, "./second/ace_model.dat")




################################## Trait Recosntruction Modeling for Phylogenetic PCA ########################################################
setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis")
library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("RPANDA")
library(ggplot2)
library(phylosignal)
library(factoextra)

source("Code/PLumage.ancestral.recon.R")	
source("Code/visualization.R")
WPTSC <- read.csv("WPTCS.csv", header = TRUE)
male_WPTSC <- read.csv("MT_Tree_Male_WPTCS.csv", header = TRUE)
female_WPTSC <- read.csv("MT_Tree_Female_WPTCS.csv", header = TRUE)

tree <- read.nexus("Cardinalidae.nexus")

acr_m_PC1_BM <- ace(phy=tree,x=as.numeric(as.character(male_score$PC1)),type="continuous",method="REML",model="BM")
acr_m_PC1_OU <- ace(phy=tree,x=as.numeric(as.character(male_score$PC1)),type="continuous",method="REML",model="OU")
acr_m_PC1_EB <- ace(phy=tree,x=as.numeric(as.character(male_score$PC1)),type="continuous",method="REML",model="EB")


#### Extract columns,
f_PC1 <- All[, "Female_pPC1"]
f_PC2 <- All[, "Female_pPC2"]
m_PC1 <- All[, "Male_pPC1"]
m_PC2 <- All[, "Male_pPC2"]
Habitat <- All[, "Habitat"]
Strata <- All[, "Strata"]
Mig <- All[,"Migration"]
#Raw Traits
f_brill <- female_WPTSC[, "AvgBrill"]
f_hue <- female_WPTSC[, "AvgHueDisp"]
m_brill <- male_WPTSC[, "AvgBrill"]
m_hue <- male_WPTSC[, "AvgHueDisp"]

#give them names
names(Mig) <- names(Strata) <- names(Habitat) <- names(m_PC2) <- names(m_PC1) <- names(m_brill) <- names(m_hue) <- rownames(All)
names(f_PC2) <- names(f_PC1) <- names(f_brill) <- names(f_hue) <- rownames(All)



# Male PC1 and raw traits
M1_BM = fitContinuous(Cardinal.tree, m_PC1 ,model="BM") 
M1_OU = fitContinuous(Cardinal.tree, m_PC1,model="OU") 
M1_EB = fitContinuous(Cardinal.tree, m_PC1,model="EB") 
M1_lamda = fitContinuous(Cardinal.tree, m_PC1,model="lambda") 
M1_kappa = fitContinuous(Cardinal.tree, m_PC1,model="kappa") 

M2_BM = fitContinuous(Cardinal.tree, m_PC2 ,model="BM") 
M2_OU = fitContinuous(Cardinal.tree, m_PC2,model="OU") 
M2_EB = fitContinuous(Cardinal.tree, m_PC2,model="EB") 
M2_lamda = fitContinuous(Cardinal.tree, m_PC2,model="lambda") 
M2_kappa = fitContinuous(Cardinal.tree, m_PC2,model="kappa") 


Mhue_BM = fitContinuous(Cardinal.tree, m_hue ,model="BM") 
Mhue_OU = fitContinuous(Cardinal.tree, m_hue,model="OU") 
Mhue_EB = fitContinuous(Cardinal.tree, m_hue,model="EB") 
Mhue_lamda = fitContinuous(Cardinal.tree, m_hue,model="lambda") 
Mhue_kappa = fitContinuous(Cardinal.tree, m_hue,model="kappa") 

mPC1_OU <- M1_OU$opt$aic
mPC1_BM <- M1_BM$opt$aic
mPC1_EB <- M1_EB$opt$aic
mPC1_lam <- M1_lamda$opt$aic
mPC1_kappa <- M1_kappa$opt$aic

mPC2_OU <- M2_OU$opt$aic
mPC2_BM <- M2_BM$opt$aic
mPC2_EB <- M2_EB$opt$aic
mPC2_lam <- M2_lamda$opt$aic
mPC2_kappa <- M2_kappa$opt$aic

mhue_OU <- Mhue_OU$opt$aic
mhue_BM <- Mhue_BM$opt$aic
mhue_EB <- Mhue_EB$opt$aic
mhue_lam <- Mhue_lamda$opt$aic
mhue_kappa <- Mhue_kappa$opt$aic

# Female PC1 and raw traits
f1_BM = fitContinuous(Cardinal.tree, f_PC1,model="BM") 
f1_OU = fitContinuous(Cardinal.tree, f_PC1,model="OU") 
f1_EB = fitContinuous(Cardinal.tree, f_PC1,model="EB") 
f1_lamda = fitContinuous(Cardinal.tree, f_PC1,model="lambda") 
f1_kappa = fitContinuous(Cardinal.tree, f_PC1,model="kappa") 

f2_BM = fitContinuous(Cardinal.tree, f_PC2,model="BM") 
f2_OU = fitContinuous(Cardinal.tree, f_PC2,model="OU") 
f2_EB = fitContinuous(Cardinal.tree, f_PC2,model="EB") 
f2_lamda = fitContinuous(Cardinal.tree, f_PC2,model="lambda") 
f2_kappa = fitContinuous(Cardinal.tree, f_PC2,model="kappa") 

fhue_BM = fitContinuous(Cardinal.tree, f_hue ,model="BM") 
fhue_OU = fitContinuous(Cardinal.tree, f_hue,model="OU") 
fhue_EB = fitContinuous(Cardinal.tree, f_hue,model="EB") 
fhue_lamda = fitContinuous(Cardinal.tree, f_hue,model="lambda") 
fhue_kappa = fitContinuous(Cardinal.tree, f_hue,model="kappa") 

fPC1_OU <- f1_OU$opt$aic
fPC1_BM <- f1_BM$opt$aic
fPC1_EB <- f1_EB$opt$aic
fPC1_lam <- f1_lamda$opt$aic
fPC1_kappa <- f1_kappa$opt$aic

fPC2_OU <- f2_OU$opt$aic
fPC2_BM <- f2_BM$opt$aic
fPC2_EB <- f2_EB$opt$aic
fPC2_lam <- f2_lamda$opt$aic
fPC2_kappa <- f2_kappa$opt$aic

fhue_OU <- fhue_OU$opt$aic
fhue_BM <- fhue_BM$opt$aic
fhue_EB <- fhue_EB$opt$aic
fhue_lam <- fhue_lamda$opt$aic
fhue_kappa <- fhue_kappa$opt$aic

AIC <- cbind(mPC1_BM,mPC1_EB,mPC1_lam,mPC1_OU,mPC1_kappa,mPC2_BM,mPC2_EB,mPC2_lam,mPC2_OU,mPC2_kappa,mhue_BM,mhue_lam,mhue_EB,mhue_OU,mhue_kappa,fPC1_BM,fPC1_EB,fPC1_lam,fPC1_OU,fPC1_kappa,fPC2_BM,fPC2_EB,fPC2_lam,fPC2_OU,fPC2_kappa,fhue_BM,fhue_lam,fhue_EB,fhue_OU,fhue_kappa)
AIC <-data.frame(AIC)
write.csv(AIC, "AIC.csv")




