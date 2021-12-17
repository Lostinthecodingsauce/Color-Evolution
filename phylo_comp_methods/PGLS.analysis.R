#New PGLS analysis


#### Reduced Patches Script 
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


setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis/Reduced_patch_analysis")
r_male_WPTSC <- read.csv("male_reduced_WPTCS.csv", header = TRUE)
r_female_WPTSC <- read.csv("female_reduced_WPTCS.csv", header = TRUE)
FullCardinal.tree <- read.nexus("Cardinalidae.nexus")
All <- read.csv("All_reduced.csv")
Mass <- read.csv("Bird.Mass.Strata.csv", header = TRUE)
tip <- c("Cyanocompsa_parellina")
Cardinal.tree <- drop.tip(FullCardinal.tree,tip) 

rownames(r_male_WPTSC) <- r_male_WPTSC[,1,]
r_male_WPTSC<- r_male_WPTSC[,c(2,3,4,5,6,7,8,9,10,11)]
r_male_WPTSC<-r_male_WPTSC[Cardinal.tree$tip.label, ]  
rownames(r_female_WPTSC) <- r_female_WPTSC[,1,]
r_female_WPTSC<- r_female_WPTSC[,c(2,3,4,5,6,7,8,9,10,11)]
r_female_WPTSC<-r_female_WPTSC[Cardinal.tree$tip.label, ]

##### Run pPCA #######
r_log.m_WPTSC <- log(r_male_WPTSC[,1:10])
#log transform data
r_log.f_WPTSC <- log(r_female_WPTSC[,1:10])
###################################################################################################################
fBrill <- r_log.f_WPTSC[,"AvgBrill"]
fChroma <- r_log.f_WPTSC[,"AvgChroma"]
fSpan <- r_log.f_WPTSC[,"AvgSpan"]
fVol <- r_log.f_WPTSC[,"Volume"]
mBrill <- r_log.m_WPTSC[,"AvgBrill"]
mChroma <- r_log.m_WPTSC[,"AvgChroma"]
mSpan <- r_log.m_WPTSC[,"AvgSpan"]
mVol <- r_log.m_WPTSC[,"Volume"]
mAhue <- r_log.m_WPTSC[,"AvgHueDisp"]
mMhue <- r_log.m_WPTSC[,"MaxHueDisp"]
fAhue <- r_log.f_WPTSC[,"AvgHueDisp"]
fMhue <- r_log.f_WPTSC[,"MaxHueDisp"]
names(fBrill) <- names(fChroma) <- names(fSpan) <- names(fVol) <- names(fAhue) <- names(fMhue) <- rownames(r_female_WPTSC)
names(mBrill) <- names(mChroma) <- names(mSpan)<-names(mVol)<- names(mAhue) <- names(mMhue) <- rownames(r_male_WPTSC)
### PGLS carper ##################
All <- read.csv("All_reduced.csv", header = TRUE)
All<- All[,c(1,2,3,4,5,6,7,8,9)]
All <- All[-c(11),]
All <- cbind(All,fBrill,mBrill)
#All<-All[Cardinal.tree$tip.label, ]  
library(caper)
comp.data<-comparative.data(Cardinal.tree, All, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
### Tested AICc for Kappa, delta, and lambda. Lambda was the best

a<-pgls(mPC1~fPC1, data=comp.data,lambda = "ML") ### Best
summary(a)
B<-pgls(mPC2~fPC2, data=comp.data,lambda = "ML")## Best
summary(B)
###########################################################################
###########################################################################
##### Effect of Strata
m1S<-pgls(mPC1~Strata, data=comp.data,lambda =  "ML")## Lambda best, Understory highyl sig
m2S<-pgls(mPC2~Strata, data=comp.data,kappa =  "ML")### Kappa best , Mid and open highly significant 
A<-pgls(mBrill~Strata, data=comp.data,lambda =  "ML")## Lambda, Only Understory **
#f1S<-pgls(f_PC1~Strata, data=comp.data,lambda =  "ML")
f1S<-pgls(fPC1~Strata, data=comp.data,delta =  "ML")### Lambda Best, None significant 
f2S<-pgls(fPC2~Strata, data=comp.data,lambda =  "ML")### lambda best , C**,O**,U*
fBS<-pgls(fBrill~Strata, data=comp.data,lambda =  "ML")## Lambda, Only Understory ***
summary(m1S)

do.call(cbind, lapply(m1S, summary))


A<-pgls(mPC2~Migration, data=comp.data,lambda =  "ML")##
B<-pgls(mPC2~Migration, data=comp.data,kappa =  "ML")
C<-pgls(mPC2~Migration, data=comp.data,delta =  "ML")
woof <- cbind(A$aicc,B$aicc,C$aicc)
woof
summary(A)

A<-pgls(mPC2~Strata, data=comp.data,lambda =  "ML")##
B<-pgls(mPC2~Strata, data=comp.data,kappa =  "ML")
C<-pgls(mPC2~Strata, data=comp.data,delta =  "ML")

m1H<-pgls(m_PC1~Habitat, data=comp.data,lambda =  "ML") # Best, not sig 
m2H<-pgls(m_PC2~Habitat, data=comp.data,lambda =  "ML")# not sig
f1H<-pgls(f_PC1~Habitat, data=comp.data,lambda =  "ML") # not sig 
f2H<-pgls(f_PC2~Habitat, data=comp.data,lambda =  "ML")# not sig

m1H<-pgls(m_PC1~Migration, data=comp.data,lambda =  "ML") # Best, not sig 
m2H<-pgls(m_PC2~Migration, data=comp.data,lambda =  "ML")# not sig
f1H<-pgls(f_PC1~Migration, data=comp.data,lambda =  "ML") # not sig 
f2H<-pgls(f_PC2~Migration, data=comp.data,lambda =  "ML")# not sig

m1F<-pgls(m_PC1~Forest, data=comp.data,lambda =  "ML") # Best, not sig 
m2F<-pgls(m_PC2~Forest, data=comp.data,lambda =  "ML")# Medium**, significant 
f1F<-pgls(f_PC1~Forest, data=comp.data,lambda =  "ML") # not sig 
f2F<-pgls(f_PC2~Forest, data=comp.data,lambda =  "ML")# not sig

m1HS<-pgls(m_PC1~Strata+Forest, data=comp.data,lambda =  "ML") # Best, not sig 
m2HS<-pgls(m_PC2~Strata+Forest, data=comp.data,lambda =  "ML")# Medium**, significant 
f1HS<-pgls(f_PC1~Strata+Forest, data=comp.data,lambda =  "ML") # not sig 
f2HS<-pgls(f_PC2~Strata+Forest, data=comp.data,lambda =  "ML")# not sig

summary(m1HS)
summary(m2HS)
summary(f1HS)
summary(f2HS)

### Comparing Best model 
A<-pgls(mPC1~Habitat, data=comp.data,lambda =  "ML")
B<-pgls(mPC1~Strata, data=comp.data,lambda =  "ML")
C<-pgls(mPC1~Migration, data=comp.data,lambda =  "ML")
D<-pgls(mPC1~Forest.Dependency, data=comp.data,lambda =  "ML")
E<-pgls(mPC1~Strata+Habitat, data=comp.data,lambda =  "ML")
E1<-pgls(mPC1~Strata+Migration, data=comp.data,lambda =  "ML")
G<-pgls(mPC1~Strata+Forest.Dependency, data=comp.data,lambda =  "ML")
H<-pgls(mPC1~Forest.Dependency+Migration, data=comp.data,lambda =  "ML")
I<-pgls(mPC1~Habitat+Migration, data=comp.data,lambda =  "ML")
J<-pgls(mPC1~Strata+Forest.Dependency+Migration, data=comp.data,lambda =  "ML")
K<-pgls(mPC1~Strata+Habitat+Migration, data=comp.data,lambda =  "ML")

woof <- cbind(A$aicc,B$aicc,C$aicc,D$aicc,E$aicc,E1$aicc,G$aicc,H$aicc,I$aicc,J$aicc,K$aicc)
woof
summary(A)
##### Raw traits 