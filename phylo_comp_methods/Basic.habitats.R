### Basic habitat analysis ###

### USe PGLS framework and run different models using habitat and foreging strata as major effects 
### Use male PC1 and Hue, Female PC1 and Hue,

#### female PC1 includes hue but not brilance, so run again using
### Male PC1, Brillance, -- Female PC1, brillance 

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
FullCardinal.tree <- read.nexus("Cardinalidae.nexus")
All <- read.csv("Species.Habs.Scores.csv")

#Make species rownames and remove extra cloumn
rownames(All) <- All[,1,]
All<- All[,c(2,3,4,5,6,7,8,9,10,11,12,13,14)]
## Subset males and females 
y <- c("Migration","Habitat","Strata","Male_pPC1","Male_pPC2","MaleAvgHueDisp","MaleAvgBrill","SDmean","SDmax")
x <- c("Migration","Habitat","Strata","Female_pPC1","Female_pPC2","FemaleAvgHueDisp","FemaleAvgBrill","SDmean","SDmax")
male_all <- All[y]
female_all <-All[x]




#### Males

##### This species is an outlier with a PC1 score of 6.5 . Should I continue to keep out?
tip <- c("Cyanocompsa_cyanoides")
Cardinal.tree <- drop.tip(FullCardinal.tree,tip) 




#### Extract columns,
f_PC1 <- All[, "Female_pPC1"]
f_PC2 <- All[, "Female_pPC2"]
m_PC1 <- All[, "Male_pPC1"]
m_PC2 <- All[, "Male_pPC2"]
Habitat <- All[, "Habitat"]
Strata <- All[, "Strata"]
Mig <- All[,"Migration"]
SDmean <- All[,"SDmean"]
SDmax <- All[,"SDmax"]
#Raw Traits
fbrill <- female_WPTSC[, "AvgBrill"]
f_hue <- female_WPTSC[, "AvgHueDisp"]
m_brill <- male_WPTSC[, "AvgBrill"]
m_hue <- male_WPTSC[, "AvgHueDisp"]

#give them names
names(SDmean) <- names(SDmax) <- names(Mig) <- names(Strata) <- names(Habitat) <- names(m_PC2) <- names(m_PC1) <- names(m_brill) <- names(m_hue) <- rownames(All)
names(f_PC2) <- names(f_PC1) <- names(f_brill) <- names(f_hue) <- rownames(All)



### Run different models of Evolution Under PGLS framework. 

#### MOdels to run
# Male complexity vs female complexity
# Male complexity vs habitat
# Male complexity vs rrphylo rates
# Male complexity vs 

### Use habitat and strata as different effects within framework. Fit female and male models seperatley using this methods
## 1) Habitat as main effect, strata as main effect, 


bm<-corBrownian(1, Cardinal.tree) # Random walk 
Pagel<-corPagel(1, Cardinal.tree) # Allows BM to devate slightly
OU<-corMartins(1, Cardinal.tree) # Evolving under stabilizing selection, or rubber brand
EB <- corBlomberg(1, Cardinal.tree) # (EB) Assume BM but rates accelerate (g<1) or decelerate (g>1) through time. Evolve Under adaptive Radiation? 

# Make a model, this is still comparing if male traits are evolving together. Spurufulous. Try again. Need SD values? 
# First argument is responce variable, second is predictor variable => "resp ~ pred"
##Males vs Females
PC1_BM <-gls(m_PC1~f_PC1, data=All, correlation=bm, method = "ML")
## Lamda produces the best model
PC1_lamda <-gls(m_PC1~f_PC1, data=All, correlation=Pagel, method = "ML")
PC1_OU <- gls(m_PC1~f_PC1, data=All, correlation=OU, method = "ML")
PC1_EB <- gls(m_PC1~f_PC1, data=All, correlation=corBlomberg(1, Cardinal.tree, fixed = TRUE), method = "ML")
summary(PC1_BM)
summary(PC1_lamda)
summary(PC1_OU)
summary(PC1_EB)

## Male and female PC1 and PC2 scores are not significantly correlated with one another


##### PC2 males vs females
PC2_BM <-gls(m_PC2~f_PC2, data=All, correlation=bm, method = "ML")
## Lamda produces the best model
PC2_lamda <-gls(m_PC2~f_PC2, data=All, correlation=Pagel, method = "ML")
PC2_OU <- gls(m_PC2~f_PC2, data=All, correlation=OU, method = "ML")
PC2_EB <- gls(m_PC2~f_PC2, data=All, correlation=corBlomberg(1, Cardinal.tree, fixed = TRUE), method = "ML")
summary(PC2_BM)
summary(PC2_lamda)
summary(PC2_OU)
summary(PC2_EB)





######################################## Sexual Dichromatism , males#########################################

mPC1_BM <-gls(m_PC1~SDmax, data=All, correlation=bm, method = "ML")
##  Lamda produces the best model
mPC1_lamda <-gls(m_PC1~SDmax, data=All, correlation=Pagel, method = "ML")
mPC1_OU <- gls(m_PC1~SDmax, data=All, correlation=OU, method = "ML")
mPC1_EB <- gls(m_PC1~SDmax, data=All, correlation=corBlomberg(1, Cardinal.tree, fixed = TRUE), method = "ML")
mPC1BMSD<-summary(mPC1_BM)
mPC1LamSD<-summary(mPC1_lamda)
mPC1OUSD<-summary(mPC1_OU)
mPC1EBSD<-summary(mPC1_EB)



#### Sexual Dichromatism , females
fPC1_BM <-gls(f_PC1~SDmax, data=All, correlation=bm, method = "ML")
##  produces the best model
fPC1_lamda <-gls(f_PC1~SDmax, data=All, correlation=Pagel, method = "ML")
fPC1_OU <- gls(f_PC1~SDmax, data=All, correlation=OU, method = "ML")
fPC1_EB <- gls(f_PC1~SDmax, data=All, correlation=corBlomberg(1, Cardinal.tree, fixed = TRUE), method = "ML")
fPC1BMsum<-summary(fPC1_BM)
fPC1Lamsum<-summary(fPC1_lamda)
fPC1OUsum<-summary(fPC1_OU)
fPC1EBsum<-summary(fPC1_EB)



####### SD, PC2
mPC2_BM <-gls(m_PC2~SDmax, data=All, correlation=bm, method = "ML")
##  Lamda produces the best model
mPC2_lamda <-gls(m_PC2~SDmax, data=All, correlation=Pagel, method = "ML")
mPC2_OU <- gls(m_PC2~SDmax, data=All, correlation=OU, method = "ML")
mPC2_EB <- gls(m_PC2~SDmax, data=All, correlation=corBlomberg(1, Cardinal.tree, fixed = TRUE), method = "ML")
mPC2BMSD<-summary(mPC2_BM)
mPC2LamSD<-summary(mPC2_lamda)
mPC2OUSD<-summary(mPC2_OU)
mPC2EBSD<-summary(mPC2_EB)

fPC2_BM <-gls(f_PC2~SDmax, data=All, correlation=bm, method = "ML")
##  produces the best model
fPC2_lamda <-gls(f_PC2~SDmax, data=All, correlation=Pagel, method = "ML")
fPC2_OU <- gls(f_PC2~SDmax, data=All, correlation=OU, method = "ML")
fPC2_EB <- gls(f_PC2~SDmax, data=All, correlation=corBlomberg(1, Cardinal.tree, fixed = TRUE), method = "ML")
fPC2BMsum<-summary(fPC2_BM)
fPC2Lamsum<-summary(fPC2_lamda)
fPC2OUsum<-summary(fPC2_OU)
fPC2EBsum<-summary(fPC2_EB)

fbrill_BM <-gls(f_PC2~SDmax, data=All, correlation=bm, method = "ML")
##  produces the best model
fPC2_lamda <-gls(f_PC2~SDmax, data=All, correlation=Pagel, method = "ML")
fPC2_OU <- gls(f_PC2~SDmax, data=All, correlation=OU, method = "ML")
fPC2_EB <- gls(f_PC2~SDmax, data=All, correlation=corBlomberg(1, Cardinal.tree, fixed = TRUE), method = "ML")
fPC2BMsum<-summary(fPC2_BM)
fPC2Lamsum<-summary(fPC2_lamda)
fPC2OUsum<-summary(fPC2_OU)
fPC2EBsum<-summary(fPC2_EB)





##### Habitat models with SD and habiat(h), strata (s), migration(m)
mPC1_H <-gls(m_PC1~SDmax+Habitat, data=All, correlation=Pagel, method = "ML")
mPC1_M <-gls(m_PC1~SDmax+Mig, data=All, correlation=Pagel, method = "ML")
mPC1_S <-gls(m_PC1~SDmax+Strata, data=All, correlation=Pagel, method = "ML")
mPC1_HM <-gls(m_PC1~SDmax+Habitat+Mig, data=All, correlation=Pagel, method = "ML")
mPC1_MH <-gls(m_PC1~SDmax+Mig+Habitat, data=All, correlation=Pagel, method = "ML")
mPC1_HS <-gls(m_PC1~SDmax+Habitat+Strata, data=All, correlation=Pagel, method = "ML")
mPC1_MS <-gls(m_PC1~SDmax+Mig+Strata, data=All, correlation=Pagel, method = "ML")

fPC1_H <-gls(f_PC1~SDmax+Habitat, data=All, correlation=Pagel, method = "ML")
fPC1_M <-gls(f_PC1~SDmax+Mig, data=All, correlation=Pagel, method = "ML")
fPC1_S <-gls(f_PC1~SDmax+Strata, data=All, correlation=Pagel, method = "ML")
fPC1_HM <-gls(f_PC1~SDmax+Habitat+Mig, data=All, correlation=Pagel, method = "ML")
fPC1_MH <-gls(f_PC1~SDmax+Mig+Habitat, data=All, correlation=Pagel, method = "ML")
fPC1_HS <-gls(f_PC1~SDmax+Habitat+Strata, data=All, correlation=Pagel, method = "ML")
fPC1_MS <-gls(f_PC1~SDmax+Mig+Strata, data=All, correlation=Pagel, method = "ML")


m1<-summary(mPC1_H)
m2<-summary(mPC1_M)
#best score
m3<-summary(mPC1_S)
m4<-summary(mPC1_HM)
m5<-summary(mPC1_MH)
#2nd best score
m6<-summary(mPC1_HS)
m7<-summary(mPC1_MS)


f1<-summary(fPC1_H)
#2nd best
f2<-summary(fPC1_M)
f3<-summary(fPC1_S)
f4<-summary(fPC1_HM)
f5<-summary(fPC1_MH)
f6<-summary(fPC1_HS)
#Best model
f7<-summary(fPC1_MS)




AIC <- cbind(m1,m2,m3,m4,m5,m6,m7,f1,f2,f3,f4,f5,f6,f7)
AIC <-data.frame(AIC)
write.table(AIC, file = "AIC_SD")
write.csv(AIC, "PGLS_results.csv", row.names=FALSE)




##### Habitat Males ######
### Strata has a significant effect in all models, habitat and migration had none
M_M <-gls(m_PC1~Mig, data=All, correlation=Pagel)
M_S <-gls(m_PC1~Strata, data=All, correlation=Pagel)
M_H <-gls(m_PC1~Habitat, data=All, correlation=Pagel)
M_HS <-gls(m_PC1~Habitat+Strata, data=All, correlation=Pagel)
M_HM <-gls(m_PC1~Habitat+Mig, data=All, correlation=Pagel)
M_MS <-gls(m_PC1~Mig+Strata, data=All, correlation=Pagel)
M_MH <-gls(m_PC1~Mig+Habitat, data=All, correlation=Pagel)


summary(M_M)
summary(M_S)
summary(M_H)
summary(M_HS)
summary(M_HM)
summary(M_MS)
summary(M_MH)




####### Females Hab,mig,Strata


f_M <-gls(f_PC1~Mig, data=All, correlation=Pagel)
f_S <-gls(f_PC1~Strata, data=All, correlation=Pagel)
f_H <-gls(f_PC1~Habitat, data=All, correlation=Pagel)
f_HS <-gls(f_PC1~Habitat+Strata, data=All, correlation=Pagel)
f_HM <-gls(f_PC1~Habitat+Mig, data=All, correlation=Pagel)
f_MS <-gls(f_PC1~Mig+Strata, data=All, correlation=Pagel)
f_MH <-gls(f_PC1~Mig+Habitat, data=All, correlation=Pagel)


summary(f_M)#sig
summary(f_S)#not sig
summary(f_H)#not sig
summary(f_HS)#not sig
summary(f_HM)#Hab slightly sig
summary(f_MS)#sig
summary(f_MH)

plot(f_PC1 ~ m_PC1)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])


M_hue_BM <-gls(M_hue~Habitat+Strata, data=All, correlation=bm)
M_hue_Lamda <-gls(m_hue~Habitat+Strata, data=All, correlation=Pagel)
M_hue_OU <- gls(m_hue~Habitat+Strata, data=All, correlation=OU)

summary(M_hue_BM)
summary(M_hue_lamda)
summary(M_hue_OU)

##### Habitat Females ######

F_HSM_BM <-gls(f_PC1~Habitat+Mig+Strata, data=All, correlation=bm)
F_HS_BM <-gls(f_PC1~Habitat+Strata, data=All, correlation=bm)
F_HSM_Lamda <-gls(f_PC1~Habitat+Mig+Strata, data=All, correlation=Pagel)
F_HSM_OU <- gls(f_PC1~Habitat+Mig+Strata, data=All, correlation=OU)

summary(F_HSM_BM)
summary(F_HS_BM)
summary(F_HSM_Lamda)
summary(F_HSM_OU)




################### Variables Effect on Male PLumage  ###################################
##### Habitat models with SD and habiat(h), strata (s), migration(m)
mPC1_H <-gls(m_PC1~Habitat, data=All, correlation=Pagel, method = "ML")
mPC1_M <-gls(m_PC1~Mig, data=All, correlation=Pagel, method = "ML")
mPC1_S <-gls(m_PC1~Strata, data=All, correlation=Pagel, method = "ML")
mPC1_HM <-gls(m_PC1~Habitat+Mig, data=All, correlation=Pagel, method = "ML")
mPC1_MH <-gls(m_PC1~Mig+Habitat, data=All, correlation=Pagel, method = "ML")
mPC1_HMS <-gls(m_PC1~Habitat+Mig+Strata, data=All, correlation=Pagel, method = "ML")

################### Variables Effect on Female PLumage  ###################################
fPC1_H <-gls(f_PC1~Habitat, data=All, correlation=Pagel, method = "ML")
fPC1_M <-gls(f_PC1~Mig, data=All, correlation=Pagel, method = "ML")
fPC1_S <-gls(f_PC1~Strata, data=All, correlation=Pagel, method = "ML")
fPC1_HM <-gls(f_PC1~Habitat+Mig, data=All, correlation=Pagel, method = "ML")
fPC1_MH <-gls(f_PC1~Mig+Habitat, data=All, correlation=Pagel, method = "ML")
fPC1_HMS <-gls(f_PC1~Habitat+Mig+Strata, data=All, correlation=Pagel, method = "ML")

m1<-summary(mPC1_H)
m2<-summary(mPC1_M)
m3<-summary(mPC1_S)
m4<-summary(mPC1_HM)
m5<-summary(mPC1_MH)
m6<-summary(mPC1_HMS)

f1<-summary(fPC1_H)
f2<-summary(fPC1_M)
f3<-summary(fPC1_S)
f4<-summary(fPC1_HM)
f5<-summary(fPC1_MH)
f6<-summary(fPC1_HMS)

AIC <- cbind(m1,m2,m3,m4,m5,m6,f1,f2,f3,f4,f5,f6)
AIC <-data.frame(AIC)


