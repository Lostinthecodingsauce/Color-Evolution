##### Substrate Matching PGLS##############
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
library(dplyr)
library(qwraps2)# Not using 
library(stargazer)

setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis/SignalingVS_CS")
FullCardinal.tree <- read.nexus("Cardinalidae.nexus")
#tip <- c("Cyanocompsa_parellina", "Amaurospiza_concolor","Chlorothraupis_stolzmanni","Cyanocompsa_brissonii")
tip <- c("Cyanocompsa_cyanoides","Amaurospiza_concolor")
Cardinal.tree <- drop.tip(FullCardinal.tree,tip)
what <- read.csv("All_substrateMatching.csv")

who <- read.csv("All_Amour.csv")
Card.tree <- drop.tip(FullCardinal.tree,tip)
tip <- c("Cyanocompsa_cyanoides")
#Create Data.frame
comp.data<-comparative.data(Cardinal.tree, what, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

comp.data<-comparative.data(Card.tree, who, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

##############################################################################################################
######################              Signaling Patches Analysis                     ##########################
##############################################################################################################
SigA<-pgls(mSigPC1~fSigPC1, data=comp.data,lambda = "ML") ### Best
summary(SigA)
plot(SigA)
abline(SigA)

SigB<-pgls(mSigPC2~fSigPC2, data=comp.data,lambda = "ML")## Best
summary(SigB)
plot(SigB)
### "Pheucticus_aureoventris" tested positivley as outlier 
######################################  Life History Traits ###############################################

#Migration
a1<-pgls(mSigPC1~Migration, data=comp.data,lambda = "ML") ### Not
b1<-pgls(mSigPC2~Migration, data=comp.data,lambda = "ML")## Not
c1<-pgls(fSigPC1~Migration, data=comp.data,lambda = "ML") ## yes
d1<-pgls(fSigPC2~Migration, data=comp.data,lambda = "ML")## Not
summary(a1)
summary(b1)
summary(c1)
summary(d1)

#Forest
a2<-pgls(mSigPC1~Forest.Dependency, data=comp.data,lambda = "ML") ### Not
b2<-pgls(mSigPC2~Forest.Dependency, data=comp.data,lambda = "ML")## Medium
c2<-pgls(fSigPC1~Forest.Dependency, data=comp.data,lambda = "ML") ## Not
d2<-pgls(fSigPC2~Forest.Dependency, data=comp.data,lambda = "ML")## Not
summary(a2)
summary(b2)

#Habitat
a3<-pgls(mSigPC1~Habitat, data=comp.data,lambda = "ML") ### Not  ### @nd best model selection 
b3<-pgls(mSigPC2~Habitat, data=comp.data,kappa  = "ML")## Not
c3<-pgls(fSigPC1~Habitat, data=comp.data,lambda = "ML") ## Not
d3<-pgls(fSigPC2~Habitat, data=comp.data,lambda = "ML")## Not
summary(a3)
summary(b3)
summary(c3)
summary(d3)


## Strata  
a4a<-pgls(mSigPC1~StrataA, data=comp.data,lambda = "ML")## Under** sig ## 3rd best model selection 
a4<-pgls(mSigPC1~Strata, data=comp.data,lambda = "ML")## Under* sig
b4a<-pgls(mSigPC2~StrataA, data=comp.data,lambda = "ML")## None
b4<-pgls(mSigPC2~Strata, data=comp.data,lambda = "ML")## None
c4a<-pgls(fSigPC1~StrataA, data=comp.data,lambda = "ML")## None
c4<-pgls(fSigPC1~Strata, data=comp.data,lambda = "ML")## None
d4a<-pgls(fSigPC2~StrataA, data=comp.data,lambda = "ML")## None
d4<-pgls(fSigPC2~Strata, data=comp.data,lambda = "ML")## Under*  U/C = 0.054
summary(a4a)
summary(b4a)
summary(c4a)
summary(d4a)

## Strata plus Migration
a5a<-pgls(mSigPC1~StrataA+Migration, data=comp.data,lambda = "ML")## Under sig**
a5<-pgls(mSigPC1~Strata+Migration, data=comp.data,lambda = "ML")## None
b5a<-pgls(mSigPC2~StrataA+Migration, data=comp.data,lambda = "ML")## None
b5<-pgls(mSigPC2~Strata+Migration, data=comp.data,lambda = "ML")## None
c5a<-pgls(fSigPC1~StrataA+Migration, data=comp.data,lambda = "ML")## Understory* +Migration** 
c5<-pgls(fSigPC1~Strata+Migration, data=comp.data,lambda = "ML")## None
d5a<-pgls(fSigPC2~StrataA+Migration, data=comp.data,lambda = "ML")## U/C* Migration*
d5<-pgls(fSigPC2~Strata+Migration, data=comp.data,lambda = "ML")## Understory** U/C*, Migration*
summary(a5a)
summary(a5)
## Strata plus Forest
a6a<-pgls(mSigPC1~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Under** highly sig ### Best model selection
a6<-pgls(mSigPC1~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Under highly sig
b6a<-pgls(mSigPC2~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Medium sig
b6<-pgls(mSigPC2~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Medium**, no forest*
c6a<-pgls(fSigPC1~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## not 
c6<-pgls(fSigPC1~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Not
d6a<-pgls(fSigPC2~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Not
d6<-pgls(fSigPC2~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## not


## Strata plus Habitat ==> StrataA good, Strata is singular  
a7a<-pgls(mSigPC1~StrataA+Habitat, data=comp.data,kappa = "ML")## terrestrial*, U/C, open habitat**
a7<-pgls(mSigPC1~Strata+Habitat, data=comp.data,delta = "ML")## Singular
b7a<-pgls(mSigPC2~StrataA+Habitat, data=comp.data,lambda = "ML")## None
b7<-pgls(mSigPC2~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
c7a<-pgls(fSigPC1~StrataA+Habitat, data=comp.data,lambda = "ML")## none 
c7<-pgls(fSigPC1~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
d7a<-pgls(fSigPC2~StrataA+Habitat, data=comp.data,lambda = "ML")## not
d7<-pgls(fSigPC2~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
summary(a7a)

summary(a7a)

# Habitat + Migration 
a8<-pgls(mSigPC1~Habitat+Migration, data=comp.data,lambda = "ML") ### no
b8<-pgls(mSigPC2~Habitat+Migration, data=comp.data,lambda = "ML")## no  #### Best model selection 
c8<-pgls(fSigPC1~Habitat+Migration, data=comp.data,lambda = "ML") ## no
d8<-pgls(fSigPC2~Habitat+Migration, data=comp.data,lambda = "ML")## no
summary(a8)


# Habitat + Forest 
a9<-pgls(mSigPC1~Forest.Dependency+Habitat, data=comp.data,lambda = "ML") ### no
b9<-pgls(mSigPC2~Forest.Dependency+Habitat, data=comp.data,lambda = "ML")## Medium **  ### 2nd best model selection 
c9<-pgls(fSigPC1~Forest.Dependency+Habitat, data=comp.data,lambda = "ML") ## no
d9<-pgls(fSigPC2~Forest.Dependency+Habitat, data=comp.data,lambda = "ML")## no

# Forest + Migration 
a10<-pgls(mSigPC1~Forest.Dependency+Migration, data=comp.data,lambda = "ML") ### no
b10<-pgls(mSigPC2~Forest.Dependency+Migration, data=comp.data,lambda = "ML")## Medium *
c10<-pgls(fSigPC1~Forest.Dependency+Migration, data=comp.data,lambda = "ML") ## no
d10<-pgls(fSigPC2~Forest.Dependency+Migration, data=comp.data,lambda = "ML")## no

#Forest vs nonforest
w<-pgls(mSigPC1~FUCK, data=comp.data,lambda = "ML") ### Not
x<-pgls(mSigPC2~FUCK, data=comp.data,lambda = "ML")## Medium
y<-pgls(fSigPC1~FUCK, data=comp.data,lambda = "ML") ## Not
z<-pgls(fSigPC2~FUCK, data=comp.data,lambda = "ML")## Not
summary(w)
summary(x)
summary(y)
summary(z)

#Forest vs nonforest
w1<-pgls(mSigPC1~FUCK+StrataA, data=comp.data,lambda = "ML") ### Not
x1<-pgls(mSigPC2~FUCK+StrataA, data=comp.data,lambda = "ML")## Medium
y1<-pgls(fSigPC1~FUCK+StrataA, data=comp.data,lambda = "ML") ## Not
z1<-pgls(fSigPC2~FUCK+StrataA, data=comp.data,lambda = "ML")## Not
summary(w1)
summary(x1)
summary(y1)
summary(z1)

w2<-pgls(mSigPC1~FUCK+Migration, data=comp.data,lambda = "ML") ### Not
x2<-pgls(mSigPC2~FUCK+Migration, data=comp.data,lambda = "ML")## Medium
y2<-pgls(fSigPC1~FUCK+Migration, data=comp.data,lambda = "ML") ## Not
z2<-pgls(fSigPC2~FUCK+Migration, data=comp.data,lambda = "ML")## Not
summary(w2)
summary(x2)
summary(y2)
summary(z2)

w3<-pgls(mSigPC1~FUCK+Forest.Dependency, data=comp.data,lambda = "ML") ### Not
x3<-pgls(mSigPC2~FUCK+Forest.Dependency, data=comp.data,lambda = "ML")## Medium
y3<-pgls(fSigPC1~FUCK+Forest.Dependency, data=comp.data,lambda = "ML") ## Not
z3<-pgls(fSigPC2~FUCK+Forest.Dependency, data=comp.data,lambda = "ML")## Not
summary(w3)
summary(x3)
summary(y3)
summary(z3)

w4<-pgls(mSigPC1~FUCK+Migration+StrataA, data=comp.data,lambda = "ML") ### Not
x4<-pgls(mSigPC2~FUCK+Migration+StrataA, data=comp.data,lambda = "ML")## Medium
y4<-pgls(fSigPC1~FUCK+Migration+StrataA, data=comp.data,lambda = "ML") ## Not
z4<-pgls(fSigPC2~FUCK+Migration+StrataA, data=comp.data,lambda = "ML")## Not
summary(w4)
summary(x4)
summary(y4)
summary(z4)


mSigPC1_SFM <-pgls(mSigPC1~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
mSigPC1_SHM <-pgls(mSigPC1~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(mSigPC1_SFM)
summary(mSigPC1_SHM)

mSigPC2_SFM <-pgls(mSigPC2~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
mSigPC2_SHM <-pgls(mSigPC2~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(mSigPC2_SFM)
summary(mSigPC2_SHM)

fSigPC1_SFM <-pgls(fSigPC1~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
fSigPC1_SHM <-pgls(fSigPC1~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(fSigPC1_SFM)
summary(fSigPC1_SHM)

fSigPC2_SFM <-pgls(fSigPC2~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
fSigPC2_SHM <-pgls(fSigPC2~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(fSigPC2_SFM)
summary(fSigPC2_SHM)
#1/a1)Migration, 2/a2)Forest.Dep, 3/a3) Habitat, 4/a4a)StrataA, 5/a5a)StrataA+Migration, 6/a6a)StrataA+Forest.Dep, 
# 7/a8)Habitat+Migration  8/a9) Forest.Dep+Habitat  9/a10)Forest.Dep+Migration
# 10/a4) Strata   11/a5) Strata+Migration   12/a6)Strata+Forest.Dep  
# 13/w) nonforest 14/w1) nonforest+StrataA 15/w2)nonforest+Migration  
#16/w4)nonforest+Migration+StrataA, 17/mSigPC1_SFM) StrataA+Forest.Dep+Habitat , 18) mSigPC1_SHM)trataA+Habitat+Mig 

m1aicc <- cbind(a1$aicc,a2$aicc,a3$aicc,a4a$aicc,a5a$aicc,a6a$aicc,a8$aicc,a9$aicc,a10$aicc,a4$aicc,
                a5$aicc,a6$aicc,w$aicc,w1$aicc,w2$aicc,w4$aicc,mSigPC1_SFM$aicc,mSigPC1_SHM$aicc) 

f1aicc <- cbind(c1$aicc,c2$aicc,c3$aicc,c4a$aicc,c5a$aicc,c6a$aicc,c8$aicc,c9$aicc,c10$aicc,c4$aicc,
                c5$aicc,c6$aicc,y$aicc,y1$aicc,y2$aicc,y4$aicc,fSigPC1_SFM$aicc,fSigPC1_SHM$aicc) 
### PC2
m2aicc <- cbind(b1$aicc,b2$aicc,b3$aicc,b4a$aicc,b5a$aicc,b6a$aicc,b8$aicc,b9$aicc,b10$aicc,b4$aicc,
                b5$aicc,b6$aicc,x$aicc,x1$aicc,x2$aicc,x4$aicc,mSigPC2_SFM$aicc,mSigPC2_SHM$aicc) 

f2aicc <- cbind(d1$aicc,d2$aicc,d3$aicc,d4a$aicc,d5a$aicc,d6a$aicc,d8$aicc,d9$aicc,d10$aicc,d4$aicc,
                d5$aicc,d6$aicc,z$aicc,z1$aicc,z2$aicc,z4$aicc,fSigPC2_SFM$aicc,fSigPC2_SHM$aicc) 

m1aicc ## 
m2aicc ## 
f1aicc ## 
f2aicc ##  
##############################################################################################################
######################              Counter-Shading Patches Analysis                ##########################
##############################################################################################################


CSa<-pgls(mCSPC1~fCSPC1, data=comp.data,lambda = "ML") ### Best
summary(CSa)
plot(CSa)
abline(CSa)

CSb<--pgls(fCSPC2~mCSPC2, data = comp.data, lambda =  "ML")
fCSPC2
summary(CSb)
plot(CSb)



######################################  Life History Traits ###############################################

#Migration
a1<-pgls(mCSPC1~Migration, data=comp.data,lambda = "ML") ### Not
b1<-pgls(mCSPC2~Migration, data=comp.data,lambda = "ML")## Not
c1<-pgls(fCSPC1~Migration, data=comp.data,lambda = "ML") ## yes
d1<-pgls(fCSPC2~Migration, data=comp.data,lambda = "ML")## Not
summary(a1)
summary(b1)
summary(c1)
summary(d1)

#Forest
a2<-pgls(mCSPC1~Forest.Dependency, data=comp.data,lambda = "ML") ### Not
b2<-pgls(mCSPC2~Forest.Dependency, data=comp.data,lambda = "ML")## Medium
c2<-pgls(fCSPC1~Forest.Dependency, data=comp.data,lambda = "ML") ## Not
d2<-pgls(fCSPC2~Forest.Dependency, data=comp.data,lambda = "ML")## Not
summary(a2)
summary(b2)

#Habitat
a3<-pgls(mCSPC1~Habitat, data=comp.data,lambda = "ML") ### Not  ### @nd best model selection 
b3<-pgls(mCSPC2~Habitat, data=comp.data,kappa  = "ML")## Not
c3<-pgls(fCSPC1~Habitat, data=comp.data,lambda = "ML") ## Not
d3<-pgls(fCSPC2~Habitat, data=comp.data,lambda = "ML")## Not
summary(a3)
summary(b3)
summary(c3)
summary(d3)


## Strata  
a4a<-pgls(mCSPC1~StrataA, data=comp.data,lambda = "ML")## Under** sig ## 3rd best model selection 
a4<-pgls(mCSPC1~Strata, data=comp.data,lambda = "ML")## Under* sig
b4a<-pgls(mCSPC2~StrataA, data=comp.data,lambda = "ML")## None
b4<-pgls(mCSPC2~Strata, data=comp.data,lambda = "ML")## None
c4a<-pgls(fCSPC1~StrataA, data=comp.data,lambda = "ML")## None
c4<-pgls(fCSPC1~Strata, data=comp.data,lambda = "ML")## None
d4a<-pgls(fCSPC2~StrataA, data=comp.data,lambda = "ML")## None
d4<-pgls(fCSPC2~Strata, data=comp.data,lambda = "ML")## Under*  U/C = 0.054
summary(a4a)
summary(b4a)
summary(c4a)
summary(d4a)

## Strata plus Migration
a5a<-pgls(mCSPC1~StrataA+Migration, data=comp.data,lambda = "ML")## Under sig**
a5<-pgls(mCSPC1~Strata+Migration, data=comp.data,lambda = "ML")## None
b5a<-pgls(mCSPC2~StrataA+Migration, data=comp.data,lambda = "ML")## None
b5<-pgls(mCSPC2~Strata+Migration, data=comp.data,lambda = "ML")## None
c5a<-pgls(fCSPC1~StrataA+Migration, data=comp.data,lambda = "ML")## Understory* +Migration** 
c5<-pgls(fCSPC1~Strata+Migration, data=comp.data,lambda = "ML")## None
d5a<-pgls(fCSPC2~StrataA+Migration, data=comp.data,lambda = "ML")## U/C* Migration*
d5<-pgls(fCSPC2~Strata+Migration, data=comp.data,lambda = "ML")## Understory** U/C*, Migration*
summary(a5a)
summary(a5)
## Strata plus Forest
a6a<-pgls(mCSPC1~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Under** highly sig ### Best model selection
a6<-pgls(mCSPC1~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Under highly sig
b6a<-pgls(mCSPC2~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Medium sig
b6<-pgls(mCSPC2~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Medium**, no forest*
c6a<-pgls(fCSPC1~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## not 
c6<-pgls(fCSPC1~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Not
d6a<-pgls(fCSPC2~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Not
d6<-pgls(fCSPC2~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## not


## Strata plus Habitat ==> StrataA good, Strata is singular  
a7a<-pgls(mCSPC1~StrataA+Habitat, data=comp.data,kappa = "ML")## terrestrial*, U/C, open habitat**
a7<-pgls(mCSPC1~Strata+Habitat, data=comp.data,delta = "ML")## Singular
b7a<-pgls(mCSPC2~StrataA+Habitat, data=comp.data,lambda = "ML")## None
b7<-pgls(mCSPC2~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
c7a<-pgls(fCSPC1~StrataA+Habitat, data=comp.data,lambda = "ML")## none 
c7<-pgls(fCSPC1~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
d7a<-pgls(fCSPC2~StrataA+Habitat, data=comp.data,lambda = "ML")## not
d7<-pgls(fCSPC2~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
summary(a7a)

summary(a7a)

# Habitat + Migration 
a8<-pgls(mCSPC1~Habitat+Migration, data=comp.data,lambda = "ML") ### no
b8<-pgls(mCSPC2~Habitat+Migration, data=comp.data,lambda = "ML")## no  #### Best model selection 
c8<-pgls(fCSPC1~Habitat+Migration, data=comp.data,lambda = "ML") ## no
d8<-pgls(fCSPC2~Habitat+Migration, data=comp.data,lambda = "ML")## no
summary(a8)


# Habitat + Forest 
a9<-pgls(mCSPC1~Forest.Dependency+Habitat, data=comp.data,lambda = "ML") ### no
b9<-pgls(mCSPC2~Forest.Dependency+Habitat, data=comp.data,lambda = "ML")## Medium **  ### 2nd best model selection 
c9<-pgls(fCSPC1~Forest.Dependency+Habitat, data=comp.data,lambda = "ML") ## no
d9<-pgls(fCSPC2~Forest.Dependency+Habitat, data=comp.data,lambda = "ML")## no

# Forest + Migration 
a10<-pgls(mCSPC1~Forest.Dependency+Migration, data=comp.data,lambda = "ML") ### no
b10<-pgls(mCSPC2~Forest.Dependency+Migration, data=comp.data,lambda = "ML")## Medium *
c10<-pgls(fCSPC1~Forest.Dependency+Migration, data=comp.data,lambda = "ML") ## no
d10<-pgls(fCSPC2~Forest.Dependency+Migration, data=comp.data,lambda = "ML")## no

#Forest vs nonforest
w<-pgls(mCSPC1~FUCK, data=comp.data,lambda = "ML") ### Not
x<-pgls(mCSPC2~FUCK, data=comp.data,lambda = "ML")## Medium
y<-pgls(fCSPC1~FUCK, data=comp.data,lambda = "ML") ## Not
z<-pgls(fCSPC2~FUCK, data=comp.data,lambda = "ML")## Not
summary(w)
summary(x)
summary(y)
summary(z)

#Forest vs nonforest
w1<-pgls(mCSPC1~FUCK+StrataA, data=comp.data,lambda = "ML") ### Not
x1<-pgls(mCSPC2~FUCK+StrataA, data=comp.data,lambda = "ML")## Medium
y1<-pgls(fCSPC1~FUCK+StrataA, data=comp.data,lambda = "ML") ## Not
z1<-pgls(fCSPC2~FUCK+StrataA, data=comp.data,lambda = "ML")## Not
summary(w1)
summary(x1)
summary(y1)
summary(z1)

w2<-pgls(mCSPC1~FUCK+Migration, data=comp.data,lambda = "ML") ### Not
x2<-pgls(mCSPC2~FUCK+Migration, data=comp.data,lambda = "ML")## Medium
y2<-pgls(fCSPC1~FUCK+Migration, data=comp.data,lambda = "ML") ## Not
z2<-pgls(fCSPC2~FUCK+Migration, data=comp.data,lambda = "ML")## Not
summary(w2)
summary(x2)
summary(y2)
summary(z2)

w3<-pgls(mCSPC1~FUCK+Forest.Dependency, data=comp.data,lambda = "ML") ### Not
x3<-pgls(mCSPC2~FUCK+Forest.Dependency, data=comp.data,lambda = "ML")## Medium
y3<-pgls(fCSPC1~FUCK+Forest.Dependency, data=comp.data,lambda = "ML") ## Not
z3<-pgls(fCSPC2~FUCK+Forest.Dependency, data=comp.data,lambda = "ML")## Not
summary(w3)
summary(x3)
summary(y3)
summary(z3)

w4<-pgls(mCSPC1~FUCK+Migration+StrataA, data=comp.data,lambda = "ML") ### Not
x4<-pgls(mCSPC2~FUCK+Migration+StrataA, data=comp.data,lambda = "ML")## Medium
y4<-pgls(fCSPC1~FUCK+Migration+StrataA, data=comp.data,lambda = "ML") ## Not
z4<-pgls(fCSPC2~FUCK+Migration+StrataA, data=comp.data,lambda = "ML")## Not
summary(w4)
summary(x4)
summary(y4)
summary(z4)


mCSPC1_SFM <-pgls(mCSPC1~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
mCSPC1_SHM <-pgls(mCSPC1~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(mCSPC1_SFM)
summary(mCSPC1_SHM)

mCSPC2_SFM <-pgls(mCSPC2~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
mCSPC2_SHM <-pgls(mCSPC2~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(mCSPC2_SFM)
summary(mCSPC2_SHM)

fCSPC1_SFM <-pgls(fCSPC1~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
fCSPC1_SHM <-pgls(fCSPC1~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(fCSPC1_SFM)
summary(fCSPC1_SHM)

fCSPC2_SFM <-pgls(fCSPC2~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
fCSPC2_SHM <-pgls(fCSPC2~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(fCSPC2_SFM)
summary(fCSPC2_SHM)
#1/a1)Migration, 2/a2)Forest.Dep, 3/a3) Habitat, 4/a4a)StrataA, 5/a5a)StrataA+Migration, 6/a6a)StrataA+Forest.Dep, 
# 7/a8)Habitat+Migration  8/a9) Forest.Dep+Habitat  9/a10)Forest.Dep+Migration
# 10/a4) Strata   11/a5) Strata+Migration   12/a6)Strata+Forest.Dep  
# 13/w) nonforest 14/w1) nonforest+StrataA 15/w2)nonforest+Migration  
#16/w4)nonforest+Migration+StrataA, 17/mCSPC1_SFM) StrataA+Forest.Dep+Habitat , 18) mCSPC1_SHM)trataA+Habitat+Mig 

m1aicc <- cbind(a1$aicc,a2$aicc,a3$aicc,a4a$aicc,a5a$aicc,a6a$aicc,a8$aicc,a9$aicc,a10$aicc,a4$aicc,
                a5$aicc,a6$aicc,w$aicc,w1$aicc,w2$aicc,w4$aicc,mCSPC1_SFM$aicc,mCSPC1_SHM$aicc) 

f1aicc <- cbind(c1$aicc,c2$aicc,c3$aicc,c4a$aicc,c5a$aicc,c6a$aicc,c8$aicc,c9$aicc,c10$aicc,c4$aicc,
                c5$aicc,c6$aicc,y$aicc,y1$aicc,y2$aicc,y4$aicc,fCSPC1_SFM$aicc,fCSPC1_SHM$aicc) 
### PC2
m2aicc <- cbind(b1$aicc,b2$aicc,b3$aicc,b4a$aicc,b5a$aicc,b6a$aicc,b8$aicc,b9$aicc,b10$aicc,b4$aicc,
                b5$aicc,b6$aicc,x$aicc,x1$aicc,x2$aicc,x4$aicc,mCSPC2_SFM$aicc,mCSPC2_SHM$aicc) 

f2aicc <- cbind(d1$aicc,d2$aicc,d3$aicc,d4a$aicc,d5a$aicc,d6a$aicc,d8$aicc,d9$aicc,d10$aicc,d4$aicc,
                d5$aicc,d6$aicc,z$aicc,z1$aicc,z2$aicc,z4$aicc,fCSPC2_SFM$aicc,fCSPC2_SHM$aicc) 

m1aicc ## 
m2aicc ## 
f1aicc ## 
f2aicc ##  





