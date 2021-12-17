
setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/Averaged.Spec.readings")

source("RtcsFunctions.R")	
load("bluetitss.dat")
source("visualization.R")

library(pavo)
library(phytools)
library(dplyr)
library(rgl)

Cardinals <- c("Amaur_con_aeq_F","Amaur_con_aeq_M","Amaur_con_rel_F","Amaur_con_rel_M","Amaur_moe_F","Amaur_moe_M","Card_card_F","Card_card_M","Card_pho_F","Card_pho_M","Card_sin_F","Card_sin_M","Cary_cana_fro_F","Cary_cana_fro_M","Cary_pol_F","Cary_pol_M","Chl_car_F","Chl_carr_M","Chl_oli_F","Chl_oli_M","Chl_stolz_F","Chl_stolz_M","Cyan_cyan_F","Cyan_cyan_M","Cyan_glau_F","Cyan_glau_M","Cyan_pare_F","Cyan_pare_M","Cyan_roth_F","Cyan_roth_M","Gran_pel_pel_F","Gran_pel_pel_M","Gran_sal_F","Gran_sal_M","Gran_ven_F","Gran_ven_M","Hab_atrim_F","Hab_atrim_M","Hab_cris_F","Hab_cris_M","Hab_fusc_F","Hab_fusc_M","Hab_gut_F","Hab_gut_M","Hab_rubica_F","Hab_rubica_M","Pass_amo_F","Pass_amo_M","Pass_bris_F","Pass_bris_M","Pass_caer_f","Pass_caer_m","Pass_ciris_F","Pass_ciris_M","Pass_cya_F","Pass_cya_M","Pass_lec_F","Pass_lec_M","Pass_ros_F","Pass_ros_M","Pass_veri_F","Pass_veri_M","Peri_ery_F","Peri_ery_M","Pheu_auro_F","Pheu_auro_M","Pheu_chry_aura_F","Pheu_chry_aura_M","Pheu_chry_F","Pheu_chry_M","Pheu_gas_F","Pheu_gas_M","Pheu_ludo_F","Pheu_ludo_M","Pheu_mel_F","Pheu_mel_M","Pheu_tib_F","Pheu_tib_M","Pir_bid_F","Pir_bid_M","Pir_ery_F","Pir_ery_M","Pir_flava_F","Pir_flava_M","Pir_hep_F","Pir_hep_M","Pir_leu_F","Pir_leu_M","Pir_ludo_F","Pir_ludo_M","Pir_lutea_F","Pir_lutea_M","Pir_oli_F","Pir_oli_M","Pir_rose_F","Pir_rose_M","Pir_rubra_F","Pir_rubra_M","Pir_ruceps_F","Pir_ruceps_M","Rhodo_cel_F","Rhodo_cel_M","Spiz_amer_F","Spiz_amer_M")

# Order
# Example Script
# For Loop for all plumage
# For loop for reduced patches
# For loop for signaling patches --> refs=subset(refs, select = -c(4,5,6,8,14,15,16,17,18))
# For loop for dcounter-shading patches --> refs=subset(refs, select = -c(2,3,7,9,10,11,12,13))


#Read in the reflectance spectra
refs <- read.csv("Amaur_moe_F_output.csv")
#refs <- read.csv("Pass_ciris_M_output.csv")
refs <- read.csv("Pir_ruceps_F_output.csv")
refs <- read.csv("Cyan_pare_F_output.csv")
refs <- read.csv("Pass_veri_F_output.csv")
#Get the relative stimulation values for each cone type for each patch.  It requries any reflectance spectra you would like analyzed, and ss, the spectral sensitivities you would like to use.  I have included the bluetit spectral sensitivities.  It will only use 300-700 nm, but data outside of this range can be submitted (it will delete them).  
refstims <- stim(refs,ss)

#Convert the stimulation values to Cartesian Coordinates (for later analyses)
refcart <- cartCoord(refstims)

#Converte the Cartesian Coordinates to Sphereical Coordinates (for later analyses)
refsphere <- sphereCoord(refcart)

#Calculate the maximum possible r-value for each hue (angles from the origin)
rmax <- rMax(refsphere)

#Calculate the acheived chroma for each patch (r/rmax)
acheivedr <- acheivedR(refsphere,rmax)

#Normalized brilliance for each patch (aka. brightness)
normbrill <- normBrill(refs)

#Color volume for all patches submitted (minimum convex polygon of all points in the tetracolorspace)
vol <- colorVolume(refcart)

#Hue disparity matrix (all patches compared to each other)
disp <- hueDisp(refsphere)

#Average, variance and maximum hue disparity
disp.summary <- summary.hueDisp(disp)

#Color span matrix (all patches compared to each other)
spans <- colorSpan(refcart)

#Average, variance and maximum color span
spans.summary <- summary.colorSpan(spans)

#Average chroma
avgChroma(refsphere)

#Average acheived chroma
avgAcheivedChroma(acheivedr)


#Average brilliance
avgBrill(normbrill)

#A number of summary measurements from all patches submitted (Average color span, variance in color span, maximum color span, color volume, average hue disparity, variance in hue disparity, maximum hue disparity, average brilliance, average chroma, and average acheived chroma)
summ <- summaryOfpatches(refs,ss)
spans <- spans

#Visualization, will output result outside of R
a<-plotPoints(refsphere, point.color = "black",type = "sphere", sphere.radius = 0.02)
b<-plotPoints(refsphere, point.color = "black",type = "sphere", sphere.radius = 0.02)
f_Pir_ruceps <-plotPoints(refsphere, point.color = "black",type = "sphere", sphere.radius = 0.02)
Cyan_pare_F <-plotPoints(refsphere, point.color = "black",type = "sphere", sphere.radius = 0.02)
Pass_veri_F <-plotPoints(refsphere, point.color = "black",type = "sphere", sphere.radius = 0.02)
f_Pir_ruceps
Cyan_pare_F
Pass_veri_F
##### For loop to analyze WPTCS for all species ######

WPTCS <-data.frame(matrix(nrow = 85, ncol = 10))
colnames(WPTCS) <- c("Species","AvgSpan","VarSpan","MaxSpan","Volume","AvgHueDisp","VarHueDisp","MaxHueDisp","AvgBrill","AvgChroma")

filenames <- list.files(path=getwd(),pattern="*.csv") 
numfiles <- length(filenames)  


datalist = list()

for (i in c(1:numfiles)){
  tryCatch({
    print(filenames[i])
    refs <- read.csv(filenames[i], header=TRUE)
    #Get the relative stimulation values for each cone type for each patch.  It requries any reflectance spectra you would like analyzed, and ss, the spectral sensitivities you would like to use.  I have included the bluetit spectral sensitivities.  It will only use 300-700 nm, but data outside of this range can be submitted (it will delete them).  
    refstims <- stim(refs,ss)
    
    #Convert the stimulation values to Cartesian Coordinates (for later analyses)
    refcart <- cartCoord(refstims)
    
    #Converte the Cartesian Coordinates to Sphereical Coordinates (for later analyses)
    refsphere <- sphereCoord(refcart)
    
    #Calculate the maximum possible r-value for each hue (angles from the origin)
    rmax <- rMax(refsphere)
    
    #Calculate the acheived chroma for each patch (r/rmax)
    acheivedr <- acheivedR(refsphere,rmax)
    
    #Normalized brilliance for each patch (aka. brightness)
    normbrill <- normBrill(refs)
    
    #Color volume for all patches submitted (minimum convex polygon of all points in the tetracolorspace)
    vol <- colorVolume(refcart)
    
    #Hue disparity matrix (all patches compared to each other)
    disp <- hueDisp(refsphere)
    
    #Average, variance and maximum hue disparity
    disp.summary <- summary.hueDisp(disp)
    
    #Color span matrix (all patches compared to each other)
    spans <- colorSpan(refcart)
    
    #Average, variance and maximum color span
    spans.summary <- summary.colorSpan(spans)
    
    #Average chroma
    avgChroma(refsphere)
    
    #Average acheived chroma
    avgAcheivedChroma(acheivedr)
    
    #Average brilliance
    avgBrill(normbrill)
    
    #A number of summary measurements from all patches submitted (Average color span, variance in color span, maximum color span, color volume, average hue disparity, variance in hue disparity, maximum hue disparity, average brilliance, average chroma, and average acheived chroma)
    summ <- summaryOfpatches(refs,ss)
    datalist[[i]] <- spans.summary},
    error = function(err) {
      
      print(paste("Didn't work:  ",filenames[i]))
      
    },
    warning = function(warn) {}, finally = {}
  )
}


### Comebine results and write them as csv
data <- do.call(rbind, datalist)
write.csv(data, file = "spans.csv")





##### For loop to analyze WPTCS for reduced patches ######
# Patches removed; Auriculars, Nape, Mantle, Side, Scaupulars, Greater Coverts, Secondaries 

WPTCS <-data.frame(matrix(nrow = 85, ncol = 10))
colnames(WPTCS) <- c("Species","AvgSpan","VarSpan","MaxSpan","Volume","AvgHueDisp","VarHueDisp","MaxHueDisp","AvgBrill","AvgChroma")

filenames <- list.files(path=getwd(),pattern="*.csv") 
numfiles <- length(filenames)  


datalist = list()

for (i in c(1:numfiles)){
  tryCatch({
    print(filenames[i])
    refs <- read.csv(filenames[i], header=TRUE)
    ## Remove extra measurements; include Alisons as well as belly and medium wing coverts 
    refs=subset(refs, select = -c(3,4,5,11,14,16,18))
    #Get the relative stimulation values for each cone type for each patch.  It requries any reflectance spectra you would like analyzed, and ss, the spectral sensitivities you would like to use.  I have included the bluetit spectral sensitivities.  It will only use 300-700 nm, but data outside of this range can be submitted (it will delete them).  
    refstims <- stim(refs,ss)
    
    #Convert the stimulation values to Cartesian Coordinates (for later analyses)
    refcart <- cartCoord(refstims)
    
    #Converte the Cartesian Coordinates to Sphereical Coordinates (for later analyses)
    refsphere <- sphereCoord(refcart)
    
    #Calculate the maximum possible r-value for each hue (angles from the origin)
    rmax <- rMax(refsphere)
    
    #Calculate the acheived chroma for each patch (r/rmax)
    acheivedr <- acheivedR(refsphere,rmax)
    
    #Normalized brilliance for each patch (aka. brightness)
    normbrill <- normBrill(refs)
    
    #Color volume for all patches submitted (minimum convex polygon of all points in the tetracolorspace)
    vol <- colorVolume(refcart)
    
    #Hue disparity matrix (all patches compared to each other)
    disp <- hueDisp(refsphere)
    
    #Average, variance and maximum hue disparity
    disp.summary <- summary.hueDisp(disp)
    
    #Color span matrix (all patches compared to each other)
    spans <- colorSpan(refcart)
    
    #Average, variance and maximum color span
    spans.summary <- summary.colorSpan(spans)
    
    #Average chroma
    avgChroma(refsphere)
    
    #Average acheived chroma
    avgAcheivedChroma(acheivedr)
    
    #Average brilliance
    avgBrill(normbrill)
    
    #A number of summary measurements from all patches submitted (Average color span, variance in color span, maximum color span, color volume, average hue disparity, variance in hue disparity, maximum hue disparity, average brilliance, average chroma, and average acheived chroma)
    summ <- summaryOfpatches(refs,ss)
    datalist[[i]] <- summ},
    error = function(err) {
      
      print(paste("Didn't work:  ",filenames[i]))
      
    },
    warning = function(warn) {}, finally = {}
  )
}


### Comebine results and write them as csv
data <- do.call(rbind, datalist)
write.csv(data, file = "reduced_WPTCS.csv")


####### Tricolor values


trivalues <- tristimulus(refs, ss)


#### Tricolor forloop #####

datalist = list()

for (i in c(1:numfiles)){
  tryCatch({
    print(filenames[i])
    refs <- read.csv(filenames[i], header=TRUE)
    trivalues <- tristimulus(refs, ss)
    datalist[[i]] <- trivalues},
    error = function(err) {
      
      print(paste("Didn't work:  ",filenames[i]))
      
    },
    warning = function(warn) {}, finally = {}
    )
}

data <- do.call(rbind, datalist)
write.csv(data, file = "All_tristimulous_Color_values.csv")




###### Sexual dichromatism attempt

library(pavo)
#Read in the reflectance spectra
refs <- read.csv("Amaur_moe_F_output.csv")
refs2 <- read.csv("Amaur_moe_M_output.csv")

refs <- read.csv("refs.csv")


####### I need to find away to match up particular patches with one another in massive matrix 

filenames <- list.files(path=getwd(),pattern="*.csv") 
numfiles <- length(filenames)  


data <- rbindlist(lapply(filenames,fread))
data <- subset(data, species == filenames[i])



All <- lapply(filenames,function(i){
    read.csv(i, header=FALSE, skip=4)
   })
data <- subset(All, by: filenames[i])


datalist = list()

for (i in c(1:numfiles)){
  tryCatch({
    print(filenames[i])
    refs <- read.csv(filenames[i], header=TRUE)



refstims <- stim(refs,ss)

#Convert the stimulation values to Cartesian Coordinates (for later analyses)
refcart <- cartCoord(refstims)

#Converte the Cartesian Coordinates to Sphereical Coordinates (for later analyses)
refsphere <- sphereCoord(refcart)

#Calculate the maximum possible r-value for each hue (angles from the origin)
rmax <- rMax(refsphere)

#Calculate the acheived chroma for each patch (r/rmax)
acheivedr <- acheivedR(refsphere,rmax)

#Normalized brilliance for each patch (aka. brightness)
normbrill <- normBrill(refs)

#Color volume for all patches submitted (minimum convex polygon of all points in the tetracolorspace)
vol <- colorVolume(refcart)

#Hue disparity matrix (all patches compared to each other)
disp <- hueDisp(refsphere)

#Average, variance and maximum hue disparity
disp.summary <- summary.hueDisp(disp)

#Color span matrix (all patches compared to each other)
spans <- colorSpan(refcart)

#Average, variance and maximum color span
spans.summary <- summary.colorSpan(spans)

#Average chroma
avgChroma(refsphere)

#Average acheived chroma
avgAcheivedChroma(acheivedr)

#Average brilliance
avgBrill(normbrill)

#A number of summary measurements from all patches submitted (Average color span, variance in color span, maximum color span, color volume, average hue disparity, variance in hue disparity, maximum hue disparity, average brilliance, average chroma, and average acheived chroma)
summ <- summaryOfpatches(refs,ss)
datalist[[i]] <- spans},
error = function(err) {
  
  print(paste("Didn't work:  ",filenames[i]))
  
},
warning = function(warn) {}, finally = {}
)
}
