### Individual plumage patch evolution
library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("RPANDA")
library(ggplot2)
library(phylosignal)
library(factoextra)
library(caper)
##############################
##############################
## Before continuing must remove all patches with brillance bellow .05. All patches below this (144 total) have 0 for other traits 
##############################
##############################
setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis")
Tri_stim <- read.csv("MT_Tree_All_tristimulous_Color_values.csv")
All <- read.csv("ALL.csv")
Cardinal.tree <- read.nexus("Cardinalidae.nexus")

bm<-corBrownian(1, Cardinal.tree) # Random walk 
Pagel<-corPagel(1, Cardinal.tree) # Allows BM to devate slightly
OU<-corMartins(1, Cardinal.tree) # Evolving under stabilizing selection, or rubber brand
EB <- corBlomberg(1, Cardinal.tree) # (EB) Assume BM but rates accelerate (g<1) or decelerate (g>1) through time. Evolve Under adaptive Radiation? 

### remove all patches with brillance bellow .05. All patches below this (110 total) have 0 for other traits
d<-Tri_stim[!(Tri_stim$normbrill <= "0.05"),]

# Seperate Tri-stimulous values by male and female
F_tri <- subset(d, Sex == "F")
M_tri <- subset(d, Sex == "M")
# Without removing bad patches 
F_tri <- subset(Tri_stim, Sex == "F")
M_tri <- subset(Tri_stim, Sex == "M")

#### Extract columns
Habitat <- All[, "Habitat"]
Strata <- All[, "Strata"]
Mig <- All[,"Migration"]
Forest <- All[,"Forest.Dependency"]

#give them names
names(Mig) <- names(Strata) <- names(Habitat) <- names(Forest) <- rownames(All)


### Strata
#A<-pgls(r~Strata, data=comp.data1,lambda = "ML") ### 
M_S<-pgls(normbrill~Strata, data=comp.data1,lambda = "ML") #
#C<-pgls(theta~Strata, data=comp.data1,lambda = "ML") #
#D<-pgls(phi~Strata, data=comp.data1,lambda = "ML")#
#summary(A)
summary(M_S) #not
#summary(C)
#summary(D)

#A<-pgls(r~Strata, data=comp.data2,lambda = "ML") ### 
F_S<-pgls(normbrill~Strata, data=comp.data2,lambda = "ML") 
#C<-pgls(theta~Strata, data=comp.data2,lambda = "ML") 
#D<-pgls(phi~Strata, data=comp.data2,lambda = "ML")#
#E<-pgls(Strata~phi+theta, data=comp.data2,lambda = "ML")#
#summary(A)
summary(F_S) #Understory**
#summary(C)
#summary(D)

######## Forest Dependency 
#A<-pgls(r~Forest, data=comp.data1,lambda = "ML") ### 
M_F<-pgls(normbrill~Forest, data=comp.data1,lambda = "ML") #
#C<-pgls(theta~Forest, data=comp.data1,lambda = "ML") #
#D<-pgls(phi~Forest, data=comp.data1,lambda = "ML")#
#summary(A)
summary(M_F)  
#summary(C)
#summary(D)

#A<-pgls(r~Forest, data=comp.data2,lambda = "ML") ### 
F_F<-pgls(normbrill~Forest, data=comp.data2,lambda = "ML") 
#C<-pgls(theta~Forest, data=comp.data2,lambda = "ML") 
#D<-pgls(phi~Forest, data=comp.data2,lambda = "ML")#
#summary(A)
summary(F_F) #Not 
#summary(C)
#summary(D)

##### Habitat
M_H<-pgls(normbrill~Habitat, data=comp.data1,lambda = "ML") #
F_H<-pgls(normbrill~Habitat, data=comp.data2,lambda = "ML") #
summary(M_H)#Not 
summary(F_H)#Not 
#### Signaling patches = Throat, Breat, Belly,
## Crown
# Male = 
# Female = normbrill = U***
mCrown <- subset(M_tri, Region == "Crown")
fCrown <- subset(F_tri, Region == "Crown")
mCrown <- cbind(mCrown,Mig,Habitat,Strata,Forest)
fCrown <- cbind(fCrown,Mig,Habitat,Strata,Forest)
### remove all patches with brillance bellow .05. All patches below this (144 total) have 0 for other traits
#mCrown<mCrown[!(mCrown$normbrill <= "0.05"),]
#fCrown<-fCrown[!(fCrown$normbrill <= "0.05"),]
#Create data frame for malE AND FEMALE PLUMAGE PATCH
comp.data1<-comparative.data(Cardinal.tree, mCrown, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
comp.data2<-comparative.data(Cardinal.tree, fCrown, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
#comp.data1<-comp.data1[Cardinal.tree$tip.label, ] 

### Breast 
# male = 
# female =     normbrill = U***  
mBreast <- subset(M_tri, Region == "Breast")
mBreast <- cbind(mBreast,Mig,Habitat,Strata,Forest)
fBreast <- subset(F_tri, Region == "Breast")
fBreast <- cbind(fBreast,Mig,Habitat,Strata,Forest)
comp.data1<-comparative.data(Cardinal.tree, mBreast, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
comp.data2<-comparative.data(Cardinal.tree, fBreast, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

### Throat 
#Male = r = All Strata**, normbrill = not signficant! 
# Female = r = All Strata**, normbrill = not signficant!  
mThroat <- subset(M_tri, Region == "Throat")
mThroat <- cbind(mThroat,Mig,Habitat,Strata,Forest)
fThroat <- subset(F_tri, Region == "Throat")
fThroat <- cbind(fThroat,Mig,Habitat,Strata,Forest)
comp.data1<-comparative.data(Cardinal.tree, mThroat, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
comp.data2<-comparative.data(Cardinal.tree, fThroat, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

# Counter-shading/sctructual patches  = Mantle, Back, Rump, Dorsal Trail

#bACK
# Male =NOT ASSOICATED WITH sTRATA, 
#female r*
mback <- subset(M_tri, Region == "Back")
mback <- cbind(mback,Mig,Habitat,Strata,Forest)
fback <- subset(F_tri, Region == "Back")
fback <- cbind(fback,Mig,Habitat,Strata,Forest)
comp.data1<-comparative.data(Cardinal.tree, mback, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
comp.data2<-comparative.data(Cardinal.tree, fback, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

#Rump
## mALE= r and Brillance ***** for all strata 
## FEMALE= r * for Mid
mRump <- subset(M_tri, Region == "Rump")
mRump <- cbind(mRump,Mig,Habitat,Strata,Forest)
fRump <- subset(F_tri, Region == "Rump")
fRump <- cbind(fRump,Mig,Habitat,Strata,Forest)
comp.data1<-comparative.data(Cardinal.tree, mRump, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
comp.data2<-comparative.data(Cardinal.tree, fRump, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

#Tail 
#Male= All not signficantly assoicated wiht traits 
# Female = All not signficantly assoicated wiht traits 
mDorsal.tail <- subset(M_tri, Region == "Dorsal.tail")
mDorsal.tail <- cbind(mDorsal.tail,Mig,Habitat,Strata,Forest)
fDorsal.tail <- subset(F_tri, Region == "Dorsal.tail")
fDorsal.tail <- cbind(fDorsal.tail,Mig,Habitat,Strata,Forest)
comp.data1<-comparative.data(Cardinal.tree, mDorsal.tail, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
comp.data2<-comparative.data(Cardinal.tree, fDorsal.tail, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

#Primaries
#Male= Brillance assoiated with open Strata , not assoicated with forest 
# Female = Not assoicated with strata or forest 
mPrimaries <- subset(M_tri, Region == "Primaries")
mPrimaries <- cbind(mPrimaries,Mig,Habitat,Strata,Forest)
fPrimaries <- subset(F_tri, Region == "Primaries")
fPrimaries <- cbind(fPrimaries,Mig,Habitat,Strata,Forest)
comp.data1<-comparative.data(Cardinal.tree, mPrimaries, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)
comp.data2<-comparative.data(Cardinal.tree, fPrimaries, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)







Crown <- cbind(mCrown, fCrown)


#### Comparing male and femlae patches
Crown_BM <-gls(mCrown~fCrown, data=Crown, correlation=bm, method = "ML")
Crown_lambda <-gls(mCrown~fCrown, data=Crown, correlation=Pagel, method = "ML")
Crown_OU <-gls(mCrown~fCrown, data=Crown, correlation=OU, method = "ML")
Crown_EB <-gls(mCrown~fCrown, data=Crown, correlation=EB, method = "ML")

summary(Crown_BM)
summary(Crown_lamda) # Significantly assoicated 
summary(Crown_OU)
summary(Crown_EB)


####### Life history vs single patch #########

m1_M <-gls(r~Mig, data=mCrown, correlation=Pagel)
m1_S <-gls(r~Strata, data=mCrown, correlation=Pagel)
m1_H <-gls(r~Habitat, data=mCrown, correlation=Pagel)
m1_HS <-gls(r~Habitat+Strata, data=mCrown, correlation=Pagel)
m1_HM <-gls(r~Habitat+Mig, data=mCrown, correlation=Pagel)
m1_MS <-gls(r~Mig+Strata, data=mCrown, correlation=Pagel)
m1_MH <-gls(r~Mig+Habitat, data=mCrown, correlation=Pagel)

summary(m1_M)
summary(m1_S)#
summary(m1_H)#
summary(m1_HS)# 
summary(m1_HM)#
summary(m1_MS)#
summary(m1_MH)


anova(m1_M,m1_S,m1_H,m1_HM,m1_MS,m1_MH)
AIC(m1_M,m1_S,m1_H,m1_HM,m1_MS,m1_MH)

m2_M <-gls(normbrill~Mig, data=mCrown, correlation=Pagel)
m2_S <-gls(normbrill~Strata, data=mCrown, correlation=Pagel)
m2_H <-gls(normbrill~Habitat, data=mCrown, correlation=Pagel)
m2_HS <-gls(normbrill~Habitat+Strata, data=mCrown, correlation=Pagel)
m2_HM <-gls(normbrill~Habitat+Mig, data=mCrown, correlation=Pagel)
m2_MS <-gls(normbrill~Mig+Strata, data=mCrown, correlation=Pagel)
m2_MH <-gls(normbrill~Mig+Habitat, data=mCrown, correlation=Pagel)
summary(m2_M)
summary(m2_S)#
summary(m2_H)#
summary(m2_HS)# 
summary(m2_HM)#
summary(m2_MS)#
summary(m2_MH)
AIC(m2_M,m2_S,m2_H,m2_HM,m2_MS,m2_MH)




#### Tests for ventral patches 
m2_M <-gls(normbrill~Mig, data=mBreast, correlation=Pagel)
m2_S <-gls(normbrill~Strata, data=mBreast, correlation=Pagel)
m2_H <-gls(normbrill~Habitat, data=mBreast, correlation=Pagel)
m2_HS <-gls(normbrill~Habitat+Strata, data=mBreast, correlation=Pagel)
m2_HM <-gls(normbrill~Habitat+Mig, data=mBreast, correlation=Pagel)
m2_MS <-gls(normbrill~Mig+Strata, data=mBreast, correlation=Pagel)
m2_MH <-gls(normbrill~Mig+Habitat, data=mBreast, correlation=Pagel)
summary(m2_M)
summary(m2_S)#
summary(m2_H)#
summary(m2_HS)# 
summary(m2_HM)#
summary(m2_MS)#
summary(m2_MH)
AIC(m2_M,m2_S,m2_H,m2_HM,m2_MS,m2_MH)




##### Tests for Male dorsal patches 

m2_M <-gls(normbrill~Mig, data=mback, correlation=Pagel)
m2_S <-gls(normbrill~Strata, data=mback, correlation=Pagel)
m2_H <-gls(normbrill~Habitat, data=mback, correlation=Pagel)
m2_HS <-gls(normbrill~Habitat+Strata, data=mback, correlation=Pagel)
m2_HM <-gls(normbrill~Habitat+Mig, data=mback, correlation=Pagel)
m2_MS <-gls(normbrill~Mig+Strata, data=mback, correlation=Pagel)
m2_MH <-gls(normbrill~Mig+Habitat, data=mback, correlation=Pagel)
summary(m2_M)
summary(m2_S)#
summary(m2_H)#
summary(m2_HS)# 
summary(m2_HM)#
summary(m2_MS)#
summary(m2_MH)
AIC(m2_M,m2_S,m2_H,m2_HM,m2_MS,m2_MH)

m2_M <-gls(normbrill~Mig, data=mRump, correlation=Pagel)
m2_S <-gls(normbrill~Strata, data=mRump, correlation=Pagel)
m2_H <-gls(normbrill~Habitat, data=mRump, correlation=Pagel)
m2_HS <-gls(normbrill~Habitat+Strata, data=mRump, correlation=Pagel)
m2_HM <-gls(normbrill~Habitat+Mig, data=mRump, correlation=Pagel)
m2_MS <-gls(normbrill~Mig+Strata, data=mRump, correlation=Pagel)
m2_MH <-gls(normbrill~Mig+Habitat, data=mRump, correlation=Pagel)
summary(m2_M)
summary(m2_S)#
summary(m2_H)#
summary(m2_HS)# 
summary(m2_HM)#
summary(m2_MS)#
summary(m2_MH)
AIC(m2_M,m2_S,m2_H,m2_HM,m2_MS,m2_MH)

#### for_loop for all patches 
datalist = list(M_tri)

for (i in c(1:datalist)){
  tryCatch({
    print(datalist[i])
    separate(speices_code, c("i..Species","Region"), sep = -1)})
}

    
    

