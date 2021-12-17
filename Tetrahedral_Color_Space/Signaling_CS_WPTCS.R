#### Signaling patches and counter shading patches script 
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
# For loop for signaling patches --> refs=subset(refs, select = c(1,2,3,9,10,11,12,13))
# Crown, Aurciculars, Throat, Breast, Side, Belly, Undertail coverts
# For loop for dcounter-shading patches --> refs=subset(refs, select = c(1,4,5,6,7,8,17)))
# Neck,Mantle, Back, Rump. Dorsal.Tail,  Primaries 
#Test
refs <- read.csv("Amaur_moe_F_output.csv")
refs=subset(refs, select = c(1,2,3,9,10,11,12,13))


WPTCS <-data.frame(matrix(nrow = 85, ncol = 10))
colnames(WPTCS) <- c("Species","AvgSpan","VarSpan","MaxSpan","Volume","AvgHueDisp","VarHueDisp","MaxHueDisp","AvgBrill","AvgChroma")

filenames <- list.files(path=getwd(),pattern="*.csv") 
numfiles <- length(filenames)  


datalist = list()

for (i in c(1:numfiles)){
  tryCatch({
    print(filenames[i])
    refs <- read.csv(filenames[i], header=TRUE)
    #For Signaling Patches
    refs=subset(refs, select = c(1,2,3,9,10,11,12,13))
    ### Switch out for Counter-shading patches 
    #refs=subset(refs, select = c(1,4,5,6,7,8,17))
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
write.csv(data, file = "SignalingPatch_WPTCS.csv")
write.csv(data, file = "CounterShadingPatch_WPTCS.csv")

##############################################################################################################
######################                          Phylogenetic PCA                    ##########################
##############################################################################################################

library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("RPANDA")
library("caper")
library(ggplot2)
library(phylosignal)
library(factoextra)
library("evomap")
library("ggpubr")
library(caper)

setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis/SignalingVS_CS")
tip <- c("Cyanocompsa_cyanoides")
#FullCard.tree <- read.nexus("MCC.Cardinalidae.nexus")
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- drop.tip(FullCard.tree,tip)
Sig.WPTCS <- read.csv(file = "SignalingPatch_WPTCS.csv")
CS.WPTCS <- read.csv(file = "CounterShadingPatch_WPTCS.csv")


##############################################################################################################
######################                          Signaling Patches PCA               ##########################
MSig.WPTCS <- subset(Sig.WPTCS, sex == "male")
FSig.WPTCS <- subset(Sig.WPTCS, sex == "female")
#Make species rownames and remove extra cloumn
rownames(MSig.WPTCS) <- MSig.WPTCS[,1,]
MSig.WPTCS<- MSig.WPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
MSig.WPTCS<-MSig.WPTCS[Card.tree$tip.label, ] 
rownames(FSig.WPTCS) <- FSig.WPTCS[,1,]
FSig.WPTCS<- FSig.WPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
FSig.WPTCS<-FSig.WPTCS[Card.tree$tip.label, ] 

#Log-transform 
log.fsig_WPTSC <- log(FSig.WPTCS[,1:10])
log.msig_WPTSC <- log(MSig.WPTCS[,1:10])

#### male PCA 
SmPCA <- phyl.pca(Card.tree, log.msig_WPTSC, method = "lambda", mode = "cor")
summary(SmPCA)#Importance of Components
print(SmPCA) #Standard Deviations
m_sig_score <- as.data.frame(SmPCA$S) ## Get scores
m_sig_loadings <- as.data.frame(SmPCA$L) # Loadings
#m_loadings <- m_loadings*-1

#### female PCA 
FPCA <- phyl.pca(Card.tree, log.fsig_WPTSC, method = "lambda", mode = "cor")
summary(FPCA)#Importance of Components
print(FPCA) #Standard Deviations
f_sig_score <- as.data.frame(FPCA$S)## Get scores
f_sig_loadings <- as.data.frame(FPCA$L) # Loadings
#f_loadings <- f_loadings*-1




###### #######################CounterShading PCA###########################
mCS.WPTCS <- subset(CS.WPTCS, sex == "male")
fCS.WPTCS <- subset(CS.WPTCS, sex == "female")
#Make species rownames and remove extra cloumn
rownames(mCS.WPTCS) <- mCS.WPTCS[,1,]
mCS.WPTCS<- mCS.WPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
mCS.WPTCS<-mCS.WPTCS[Card.tree$tip.label, ] 
rownames(fCS.WPTCS) <- fCS.WPTCS[,1,]
fCS.WPTCS<- fCS.WPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
fCS.WPTCS<-fCS.WPTCS[Card.tree$tip.label, ] 
# Replace any birds with strange volumes (0 or infinte) with volume of 0.000000001
mCS.WPTCS$Volume <- as.numeric(gsub(pattern = '^0$', replacement = 0.000000001, mCS.WPTCS$Volume))
fCS.WPTCS$Volume <- as.numeric(gsub(pattern = '^0$', replacement = 0.000000001, fCS.WPTCS$Volume))
#Log-transform 
log.fCS_WPTSC <- log(fCS.WPTCS[,1:10])
log.mCS_WPTSC <- log(mCS.WPTCS[,1:10])

#### male PCA 
CSmPCA <- phyl.pca(Card.tree, log.mCS_WPTSC, method = "lambda", mode = "cor")
summary(CSmPCA)#Importance of Components
print(CSmPCA) #Standard Deviations
m_CS_score <- as.data.frame(CSmPCA$S) ## Get scores
m_CS_loadings <- as.data.frame(CSmPCA$L) # Loadings
#m_loadings <- m_loadings*-1

#### female PCA 
CSFPCA <- phyl.pca(Card.tree, log.fCS_WPTSC, method = "lambda", mode = "cor")
summary(CSFPCA)#Importance of Components
print(CSFPCA) #Standard Deviations
f_CS_score <- as.data.frame(CSFPCA$S)## Get scores
f_CS_loadings <- as.data.frame(CSFPCA$L) # Loadings
#f_loadings <- f_loadings*-1



####### Plotting
All <- read.csv("All_substrateMatching.csv")
male <- ggplot(data = All, aes(x = mSigPC1, y = mSigPC2, label = rownames(All))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Signaling PC1") + geom_point(aes(size = 4)) #(colour=All$Genus,size = 4)) 
male + theme_classic() 

female <- ggplot(data = All, aes(x = fSigPC1, y = fSigPC2, label = rownames(All))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Female Sig PC1 v PC2") + geom_point(aes(size = 4)) 
female + theme_classic() 

male <- ggplot(data = All, aes(x = mCSPC1, y = mCSPC2, label = rownames(All))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Male CS PC1vPC2") + geom_point(aes(size = 4)) #(colour=All$Genus,size = 4)) 
male + theme_classic() 
female <- ggplot(data = All, aes(x = fCSPC1, y = fCSPC2, label = rownames(All))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Female CS PC1 v PC2") + geom_point(aes(size = 4)) 
female + theme_classic() 

##### PGLS 

#signaling male vs female 
library(caper)
library(ggplot2)
All <- read.csv("All_substrateMatching.csv")
comp.data<-comparative.data(Card.tree, All, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

SigPC1<-pgls(mSigPC1~fSigPC1, data=comp.data,lambda =  "ML")
summary(SigPC1)
plot(SigPC1)
## Highly Significant, but plots look bad. "Pheucticus_aureoventris" was tested positivley as an outlier
SigPC2<-pgls(mSigPC2~fSigPC2, data=comp.data,lambda =  "ML")
summary(SigPC2)
plot(SigPC2)

#### Not correlated, plots looks good! No outliers
CSPC1<-pgls(mCSPC1~fCSPC1, data=comp.data,lambda =  "ML")
summary(CSPC1)
plot(CSPC1)
## Plots look bad. Highly correlated, no outliers identified 
CSPC2<-pgls(mCSPC2~fCSPC2, data=comp.data,lambda =  "ML")
summary(CSPC2)
plot(CSPC2)


######
tip <- c("Cyanocompsa_cyanoides","Amaurospiza_concolor")
Card.tree <- drop.tip(FullCardinal.tree,tip)
Butt <- read.csv("All_substrateMatching.csv")
rownames(Butt) <- Butt[,1,]
Butt<- Butt[,c(2,3,4,5,6,7,8,9,10,11,12,13,14)]
Butt<- Butt[Card.tree$tip.label, ]

bm<-corBrownian(1, Card.tree) # Random walk 
Pagel<-corPagel(1, Card.tree) # Allows BM to devate slightly
OU<-corMartins(1, Card.tree) # Evolving under stabilizing selection, or rubber brand
EB <- corBlomberg(1, Card.tree) # (EB) Assume BM but rates accelerate (g<1) or decelerate (g>1) through time. Evolve Under adaptive Radiation? 

# PC1
M_F <-gls(mSigPC1~fSigPC1, data=Butt, correlation=Pagel)
summary(M_F) # Correlated ***
plot(mSigPC1~fSigPC1, data = Butt)
abline(M_F)

M_F <-gls(mSigPC2~fSigPC2, data=Butt, correlation=Pagel)
summary(M_F) # Correlated ***
plot(mSigPC2~fSigPC2, data = Butt)
abline(M_F)

########
######## Counter-shading
########
# PC1
M_F <-gls(mCSPC1~fCSPC1, data=Butt, correlation=Pagel)
summary(M_F) # Not correlated at all!!!!
plot(mCSPC1~fCSPC1, data = Butt)
abline(M_F)

M_F2 <-gls(mCSPC2~fCSPC2, data=Butt, correlation=Pagel)
summary(M_F2) # Correlated ****
plot(mCSPC2~fCSPC2, data = Butt)
abline(M_F2)
