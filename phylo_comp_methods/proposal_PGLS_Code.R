###### Analysis of life history traits impact on plumage evolution
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
library(BTprocessR)

setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis")
#tip <- c("Caryothraustes_canadensis","Cyanocompsa_cyanoides")
tip <- c("Cyanocompsa_cyanoides")
#FullCard.tree <- read.nexus("MCC.Cardinalidae.nexus")
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- drop.tip(FullCard.tree,tip)


#####

This <- read.csv("spatial.data.csv", header = TRUE) # without Cyanocompsa cyanoides
#Wack <- read.csv("wack.csv", header = TRUE) # without Caryothraustes_canadensis


## Build comparative data frame  

comp.data<-comparative.data(FullCard.tree, This, names.col = "ï..species", vcv.dim=2, warn.dropped=TRUE)

#comp.data<-comparative.data(FullCard.tree, Wack, names.col = "ï..species", vcv.dim=2, warn.dropped=TRUE)

###########
a<-pgls(fPC1~mPC1, data=comp.data, lambda =  "ML") ### Best
summary(a)
plot(a)
abline(a)
lambda.profile1 <- pgls.profile(a, "lambda" )
plot(lambda.profile1)## Profile is kinda flat, not quite an optimum peak. Compare to OLS (lambda = 0) and INdy contrast (Lambda=1)


B<-pgls(fPC2~mPC2, data=comp.data, lambda = "ML")## Best
summary(B)
plot(B)
lambda.profile2 <- pgls.profile(B, "lambda" )
plot(lambda.profile2)

mPC1vfPC1<- format.pval(pf(oof$fstatistic[1], oof$fstatistic[2], oof$fstatistic[3], lower.tail = FALSE))
mPC2vfPC2<- format.pval(pf(oof$fstatistic[1], oof$fstatistic[2], oof$fstatistic[3], lower.tail = FALSE))
waaf <- cbind(mPC1vfPC1,mPC2vfPC2)

#Migration
a1<-pgls(mPC1~Migration, data=comp.data,lambda = "ML") ### Not
b1<-pgls(mPC2~Migration, data=comp.data,lambda = "ML")## Not
c1<-pgls(fPC1~Migration, data=comp.data,lambda = "ML") ## yes
d1<-pgls(fPC2~Migration, data=comp.data,lambda = "ML")## Not
summary(a1)
summary(b1)
summary(c1)
summary(d1)

#Forest
a2<-pgls(mPC1A~Forest.Dependency, data=comp.data,lambda = "ML") ### Not
b2<-pgls(mPC2A~Forest.Dependency, data=comp.data,lambda = "ML")## Medium
c2<-pgls(fPC1A~Forest.Dependency, data=comp.data,lambda = "ML") ## Not
d2<-pgls(fPC2A~Forest.Dependency, data=comp.data,lambda = "ML")## Not
summary(a2)
summary(b2)

#Habitat
a3<-pgls(mPC1A~Habitat, data=comp.data,lambda = "ML") ### Not  ### @nd best model selection 
b3<-pgls(mPC2A~Habitat, data=comp.data,kappa  = "ML")## Not
c3<-pgls(fPC1A~Habitat, data=comp.data,lambda = "ML") ## Not
d3<-pgls(fPC2A~Habitat, data=comp.data,lambda = "ML")## Not
summary(a3)
summary(b3)
summary(c3)
summary(d3)


## Strata  
a4a<-pgls(mPC1A~StrataA, data=comp.data,lambda = "ML")## Under** sig ## 3rd best model selection 
a4<-pgls(mPC1A~Strata, data=comp.data,lambda = "ML")## Under* sig
b4a<-pgls(mPC2A~StrataA, data=comp.data,lambda = "ML")## None
b4<-pgls(mPC2A~Strata, data=comp.data,lambda = "ML")## None
c4a<-pgls(fPC1A~StrataA, data=comp.data,lambda = "ML")## None
c4<-pgls(fPC1A~Strata, data=comp.data,lambda = "ML")## None
d4a<-pgls(fPC2A~StrataA, data=comp.data,lambda = "ML")## None
d4<-pgls(fPC2A~Strata, data=comp.data,lambda = "ML")## Under*  U/C = 0.054
summary(a4a)
summary(b4a)
summary(c4a)
summary(d4a)

## Strata plus Migration
a5a<-pgls(mPC1A~StrataA+Migration, data=comp.data,lambda = "ML")## Under sig**
a5<-pgls(mPC1A~Strata+Migration, data=comp.data,lambda = "ML")## None
b5a<-pgls(mPC2A~StrataA+Migration, data=comp.data,lambda = "ML")## None
b5<-pgls(mPC2A~Strata+Migration, data=comp.data,lambda = "ML")## None
c5a<-pgls(fPC1A~StrataA+Migration, data=comp.data,lambda = "ML")## Understory* +Migration** 
c5<-pgls(fPC1A~Strata+Migration, data=comp.data,lambda = "ML")## None
d5a<-pgls(fPC2A~StrataA+Migration, data=comp.data,lambda = "ML")## U/C* Migration*
d5<-pgls(fPC2A~Strata+Migration, data=comp.data,lambda = "ML")## Understory** U/C*, Migration*
summary(a5a)
summary(a5)
## Strata plus Forest
a6a<-pgls(mPC1A~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Under** highly sig ### Best model selection
a6<-pgls(mPC1A~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Under highly sig
b6a<-pgls(mPC2A~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Medium sig
b6<-pgls(mPC2A~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Medium**, no forest*
c6a<-pgls(fPC1A~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## not 
c6<-pgls(fPC1A~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Not
d6a<-pgls(fPC2A~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Not
d6<-pgls(fPC2A~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## not
summary(a6a)
summary(b6a)
summary(b6a)
summary(d6a)

## Strata plus Habitat ==> StrataA good, Strata is singular  
a7a<-pgls(mPC1A~StrataA+Habitat, data=comp.data,kappa = "ML")## terrestrial*, U/C, open habitat**
a7<-pgls(mPC1A~Strata+Habitat, data=comp.data,delta = "ML")## Singular
b7a<-pgls(mPC2A~StrataA+Habitat, data=comp.data,lambda = "ML")## None
b7<-pgls(mPC2A~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
c7a<-pgls(fPC1A~StrataA+Habitat, data=comp.data,lambda = "ML")## none 
c7<-pgls(fPC1A~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
d7a<-pgls(fPC2A~StrataA+Habitat, data=comp.data,lambda = "ML")## not
d7<-pgls(fPC2A~Strata+Habitat, data=comp.data,lambda = "ML")## Singular
summary(a7a)

summary(a7a)

# Habitat + Migration 
a8<-pgls(mPC1A~Habitat+Migration, data=comp.data,lambda = "ML") ### no
b8<-pgls(mPC2A~Habitat+Migration, data=comp.data,lambda = "ML")## no  #### Best model selection 
c8<-pgls(fPC1A~Habitat+Migration, data=comp.data,lambda = "ML") ## no
d8<-pgls(fPC2A~Habitat+Migration, data=comp.data,lambda = "ML")## no
summary(c8)


# Habitat + Forest 
a9<-pgls(mPC1~Forest.Dependency+Habitat, data=comp.data,lambda = "ML") ### no
b9<-pgls(mPC2~Forest.Dependency+Habitat, data=comp.data,lambda = "ML")## Medium **  ### 2nd best model selection 
c9<-pgls(fPC1~Forest.Dependency+Habitat, data=comp.data,lambda = "ML") ## no
d9<-pgls(fPC2~Forest.Dependency+Habitat, data=comp.data,lambda = "ML")## no
summary(a9)
# Forest + Migration 
a10<-pgls(mPC1~Forest.Dependency+Migration, data=comp.data,lambda = "ML") ### no
b10<-pgls(mPC2~Forest.Dependency+Migration, data=comp.data,lambda = "ML")## Medium *
c10<-pgls(fPC1~Forest.Dependency+Migration, data=comp.data,lambda = "ML") ## no
d10<-pgls(fPC2~Forest.Dependency+Migration, data=comp.data,lambda = "ML")## no
summary(a10)
summary(b10)



#Forest vs nonforest
w<-pgls(mPC1~FUCK, data=comp.data,lambda = "ML") ### Not
x<-pgls(mPC2~FUCK, data=comp.data,lambda = "ML")## Medium
y<-pgls(fPC1~FUCK, data=comp.data,lambda = "ML") ## Not
z<-pgls(fPC2~FUCK, data=comp.data,lambda = "ML")## Not
summary(w)
summary(x)
summary(y)
summary(z)

#Forest vs nonforest
w1<-pgls(mPC1~FUCK+StrataA, data=comp.data,lambda = "ML") ### Not
x1<-pgls(mPC2~FUCK+StrataA, data=comp.data,lambda = "ML")## Medium
y1<-pgls(fPC1~FUCK+StrataA, data=comp.data,lambda = "ML") ## Not
z1<-pgls(fPC2~FUCK+StrataA, data=comp.data,lambda = "ML")## Not
summary(w1)
summary(x1)
summary(y1)
summary(z1)

w2<-pgls(mPC1~FUCK+Migration, data=comp.data,lambda = "ML") ### Not
x2<-pgls(mPC2~FUCK+Migration, data=comp.data,lambda = "ML")## Medium
y2<-pgls(fPC1~FUCK+Migration, data=comp.data,lambda = "ML") ## Not
z2<-pgls(fPC2~FUCK+Migration, data=comp.data,lambda = "ML")## Not
summary(w2)
summary(x2)
summary(y2)
summary(z2)

w3<-pgls(mPC1~FUCK+Forest.Dependency, data=comp.data,lambda = "ML") ### Not
x3<-pgls(mPC2~FUCK+Forest.Dependency, data=comp.data,lambda = "ML")## Medium
y3<-pgls(fPC1~FUCK+Forest.Dependency, data=comp.data,lambda = "ML") ## Not
z3<-pgls(fPC2~FUCK+Forest.Dependency, data=comp.data,lambda = "ML")## Not
summary(w3)
summary(x3)
summary(y3)
summary(z3)

w4<-pgls(mPC1~FUCK+Migration+StrataA, data=comp.data,lambda = "ML") ### Not
x4<-pgls(mPC2~FUCK+Migration+StrataA, data=comp.data,lambda = "ML")## Medium
y4<-pgls(fPC1~FUCK+Migration+StrataA, data=comp.data,lambda = "ML") ## Not
z4<-pgls(fPC2~FUCK+Migration+StrataA, data=comp.data,lambda = "ML")## Not
summary(w4)
summary(x4)
summary(y4)
summary(z4)

 
mPC1_SFM <-pgls(mPC1~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
mPC1_SHM <-pgls(mPC1~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(mPC1_SFM)
summary(mPC1_SHM)

mPC2_SFM <-pgls(mPC2~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
mPC2_SHM <-pgls(mPC2~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(mPC2_SFM)
summary(mPC2_SHM)

fPC1_SFM <-pgls(fPC1~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
fPC1_SHM <-pgls(fPC1~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(fPC1_SFM)
summary(fPC1_SHM)

fPC2_SFM <-pgls(fPC2~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
fPC2_SHM <-pgls(fPC2~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
summary(fPC2_SFM)
summary(fPC2_SHM)
#1/a1)Migration, 2/a2)Forest.Dep, 3/a3) Habitat, 4/a4a)StrataA, 5/a5a)StrataA+Migration, 6/a6a)StrataA+Forest.Dep, 
# 7/a8)Habitat+Migration  8/a9) Forest.Dep+Habitat  9/a10)Forest.Dep+Migration
# 10/a4) Strata   11/a5) Strata+Migration   12/a6)Strata+Forest.Dep  
# 13/w) nonforest 14/w1) nonforest+StrataA 15/w2)nonforest+Migration  
#16/w4)nonforest+Migration+StrataA, 17/mPC1_SFM) StrataA+Forest.Dep+Habitat , 18) mPC1_SHM)trataA+Habitat+Mig 

m1aicc <- cbind(a1$aicc,a2$aicc,a3$aicc,a4a$aicc,a5a$aicc,a6a$aicc,a8$aicc,a9$aicc,a10$aicc,a4$aicc,
                a5$aicc,a6$aicc,w$aicc,w1$aicc,w2$aicc,w4$aicc,mPC1_SFM$aicc,mPC1_SHM$aicc) 

f1aicc <- cbind(c1$aicc,c2$aicc,c3$aicc,c4a$aicc,c5a$aicc,c6a$aicc,c8$aicc,c9$aicc,c10$aicc,c4$aicc,
                c5$aicc,c6$aicc,y$aicc,y1$aicc,y2$aicc,y4$aicc,fPC1_SFM$aicc, fPC1_SHM$aicc) 
### PC2
m2aicc <- cbind(b1$aicc,b2$aicc,b3$aicc,b4a$aicc,b5a$aicc,b6a$aicc,b8$aicc,b9$aicc,b10$aicc,b4$aicc,
                b5$aicc,b6$aicc,x$aicc,x1$aicc,x2$aicc,x4$aicc,mPC2_SFM$aicc,mPC2_SHM$aicc) 

f2aicc <- cbind(d1$aicc,d2$aicc,d3$aicc,d4a$aicc,d5a$aicc,d6a$aicc,d8$aicc,d9$aicc,d10$aicc,d4$aicc,
                d5$aicc,d6$aicc,z$aicc,z1$aicc,z2$aicc,z4$aicc,fPC2_SFM$aicc,fPC2_SHM$aicc) 

m1aicc ## 
m2aicc ## 
f1aicc ## 
f2aicc ##  
Woof <- cbind.data.frame(m1aicc,m2aicc,f1aicc,f2aicc)
m1ANOVA<- anova(a1,a2,a3,a4a,a5a,a6a,a7a,a8,a9,a10,a4,a5,a6,mPC1_SFM,mPC1_SHM) 

#OLD data
#1/a1)Migration, 2/a2)Forest.Dep, 3/a3) Habitat, 4/a4a)StrataA, 5/a5a)StrataA+Migration, 6/a6a)StrataA+Forest.Dep, (GONE!) 7/a7a)StrataA+Habitat, 
# 8/a8)Habitat+Migration  9/a9) Forest.Dep+Habitat  10/a10)Forest.Dep+Migration
# 11/a4) Strata   12/a5) Strata+Migration   13/a6)Strata+Forest.Dep  
# 14/w) nonforest 15/w1) nonforest+StrataA 16/w2)nonforest+Migration 17/w3)nonforest+Forest.Dep ## Canot USE!
#18/w4)nonforest+Migration+StrataA, 19/mPC1_SFM) StrataA+Forest.Dep+Habitat (GONE!) 20/mPC1_SHM) StrataA+Habitat+ Migration 
f2aicc <- cbind(d1$aicc,d2$aicc,d3$aicc,d4a$aicc,d5a$aicc,d6a$aicc,d7a$aicc,d8$aicc,d9$aicc,d10$aicc,d4$aicc,
                d5$aicc,d6$aicc,z$aicc,z1$aicc,z2$aicc,z4$aicc,fPC2_SFM$aicc,fPC2_SHM$aicc) 

############################################################################################











###########################################################################
###########################################################################
##### Comparng models where males are responce vs where females are the response 
###########################################################################
###########################################################################
############ Males
M1<-pgls(MaleBrill~Migration, data=comp.data,lambda = "ML")## 
M2<-pgls(MaleBrill~Forest.Dependency, data=comp.data,lambda = "ML")## 
M3<-pgls(MaleBrill~Habitat, data=comp.data,lambda = "ML")## 
M4<-pgls(MaleBrill~StrataA, data=comp.data,lambda = "ML")## 
M5<-pgls(MaleBrill~StrataA+Migration, data=comp.data,lambda = "ML")##
M6<-pgls(MaleBrill~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")##
M7<-pgls(MaleBrill~StrataA+Habitat, data=comp.data,lambda = "ML")##
M8<-pgls(MaleBrill~Habitat+Migration, data=comp.data,lambda = "ML")##
M9<-pgls(MaleBrill~Habitat+Forest.Dependency, data=comp.data,lambda = "ML")##
M10<-pgls(MaleBrill~Forest.Dependency+Migration, data=comp.data,lambda = "ML")##
#M11<-pgls(MaleBrill~Strata, data=comp.data,lambda = "ML")##
#M12<-pgls(MaleBrill~Strata+Migration, data=comp.data,lambda = "ML")##
#M13<-pgls(MaleBrill~Strata+Forest.Dependency, data=comp.data,lambda = "ML")##
#M14<-pgls(MaleBrill~Strata+Habitat, data=comp.data,lambda = "ML")##
M11<-pgls(MaleBrill~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")##
M12<-pgls(MaleBrill~StrataA+Habitat+Forest.Dependency, data=comp.data,lambda = "ML")##
M13<-pgls(MaleBrill~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")##
M14<-pgls(MaleBrill~Habitat+Forest.Dependency+Migration, data=comp.data,lambda = "ML")##
#RESULTS
summary(M1)
summary(M2)
summary(M3)
summary(M4)# Yes under
summary(M5)# Yes under, no mig
summary(M6)
summary(M7)
summary(M8)
summary(M9)
summary(M10)
summary(M11)
summary(M12)
summary(M13)
summary(M14)
#summary(M15)
#summary(M16)
#summary(M16)
#summary(M17)
#summary(M18)
MaleBrillModels <- cbind(M1$aicc,M2$aicc,M3$aicc,M4$aicc,M5$aicc,M6$aicc,M7$aicc,M8$aicc,M9$aicc,M10$aicc,M11$aicc,M12$aicc,M13$aicc,M14$aicc)#,M15$aicc,M16$aicc,M17$aicc,M18$aicc)
############# Females ###################
############ Males
F1<-pgls(FemaleBrill~Migration, data=comp.data,lambda = "ML")## 
F2<-pgls(FemaleBrill~Forest.Dependency, data=comp.data,lambda = "ML")## 
F3<-pgls(FemaleBrill~Habitat, data=comp.data,lambda = "ML")## 
F4<-pgls(FemaleBrill~StrataA, data=comp.data,lambda = "ML")## 
F5<-pgls(FemaleBrill~StrataA+Migration, data=comp.data,lambda = "ML")##
F6<-pgls(FemaleBrill~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")##
F7<-pgls(FemaleBrill~StrataA+Habitat, data=comp.data,lambda = "ML")##
F8<-pgls(FemaleBrill~Habitat+Migration, data=comp.data,lambda = "ML")##
F9<-pgls(FemaleBrill~Habitat+Forest.Dependency, data=comp.data,lambda = "ML")##
F10<-pgls(FemaleBrill~Forest.Dependency+Migration, data=comp.data,lambda = "ML")##
#F11<-pgls(FemaleBrill~Strata, data=comp.data,lambda = "ML")##
#F12<-pgls(FemaleBrill~Strata+Migration, data=comp.data,lambda = "ML")##
#F13<-pgls(FemaleBrill~Strata+Forest.Dependency, data=comp.data,lambda = "ML")##
#F14<-pgls(FemaleBrill~Strata+Habitat, data=comp.data,lambda = "ML")##
F11<-pgls(FemaleBrill~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")##
F12<-pgls(FemaleBrill~StrataA+Habitat+Forest.Dependency, data=comp.data,lambda = "ML")##
F13<-pgls(FemaleBrill~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")##
F14<-pgls(FemaleBrill~Habitat+Forest.Dependency+Migration, data=comp.data,lambda = "ML")##

summary(F1)
summary(F2)
summary(F3)
summary(F4)#Stata under
summary(F5)#o ly stra
summary(F6)
summary(F7)
summary(F8)
summary(F9)
summary(F10)
summary(F11)
summary(F12)
summary(F13)
summary(F14)
#summary(F15)
#summary(F16)
#summary(F16)
#summary(F17)
#summary(F18)
FemaleBrillanceModels <- cbind(F1$aicc, F2$aicc,F3$aicc,F4$aicc,F5$aicc,F6$aicc,F7$aicc,F8$aicc,F9$aicc,F10$aicc,F11$aicc,F12$aicc,F13$aicc,F14$aicc) #,F15$aicc,F16$aicc,F17$aicc,F18$aicc)



########################## Other method ##################
# Extract columns
mPC1 <- Butt[, "mPC1"]
fPC1A <- Butt[, "fPC1A"]
mPC2 <- Butt[, "mPC2"]
fPC2A <- Butt[, "fPC2A"]

# Give them names
names(fPC1A) <- names(mPC1) <- rownames(Butt)

# Calculate PICs
hPic <- pic(mPC1, Card.tree)
aPic <- pic(fPC1A, Card.tree)

# Make a model
picModel <- lm(hPic ~ aPic - 1)

# Yes, significant
summary(picModel)



############
tip <- c("Cyanocompsa_cyanoides")
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Cardinal.tree <- drop.tip(FullCard.tree,tip)

#Butt <- read.csv("All.csv", header = TRUE)
#tip <- c("Cyanocompsa_cyanoides","Amaurospiza_concolor")
Butt <- read.csv("New.csv", header = TRUE)
rownames(Butt) <- Butt[,1,]
Butt<-Butt[Cardinal.tree$tip.label, ] 
Butt<- Butt[,c(2,3,4,5,6,7,8,9,10,11,12,13,14)]


#Fix error of OU with convergence
tempTree <- Cardinal.tree
tempTree$edge.length <- tempTree$edge.length * 100
##### create PGLS background 
bm<-corBrownian(1, Cardinal.tree) # Random walk 
Pagel<-corPagel(1, Cardinal.tree) # Buttows BM to devate slightly
#OU <- corMartins(1, Cardinal.tree)
OU <- corMartins(1, phy = tempTree, fixed = FALSE) # Evolving under stabilizing selection, or rubber brand
EB <- corBlomberg(1, Cardinal.tree, fixed = TRUE) # (EB) Assume BM but rates accelerate (g<1) or decelerate (g>1) through time. Evolve Under adaptive Radiation? 

### method = Maximum Likelihodd value for each ML parameter (kappa, delta, lambda) will be found within set bounds

##Males vs Females PC1. ## Lamda produces the best model
PC1_BM <-gls(mPC1~fPC1, data=Butt, correlation=bm, method = "ML")
PC1_lamda <-gls(fPC1A~mPC1, data=Butt, correlation=Pagel, method = "ML")
PC1_OU <- gls(mPC1~fPC1, data=Butt, correlation=OU, method = "ML")
PC1_EB <- gls(mPC1~fPC1, data=Butt, correlation=corBlomberg(1, Cardinal.tree, fixed = TRUE), method = "ML")
summary(PC1_lamda) # Significantly assoicated 
coef(PC1_lamda)
plot(fPC1A ~ mPC1)
abline(a = coef(PC1_lamda)[1], b = coef(PC1_lamda)[2])

anova(PC1_BM,PC1_EB,PC1_lamda,PC1_OU)
AIC(PC1_BM,PC1_EB,PC1_lamda,PC1_OU)

# PC2
PC2_lamda <-gls(fPC2A~mPC2, data=Butt, correlation=Pagel)
summary(PC2_lamda)
coef(PC2_lamda)
plot(fPC2A ~ mPC2)
abline(a = coef(PC2_lamda)[1], b = coef(PC2_lamda)[2])

###### Plotting 
M1_F1 <-gls(mPC1~fPC1A, data=Butt, correlation=Pagel,method = "ML")
plot(mPC1~fPC1A, data=Butt, correlation=Pagel,method = "ML")
abline(M1_F1)

M2_F2 <-gls(mPC2~fPC2A, data=Butt, correlation=Pagel,method = "ML")
summary(M2_F2)
plot(fPC2A~mPC2, data = comp.data,correlation=Pagel,method = "ML")
abline(M2_F2)

ggscatter(what, x = "mPC1", y = "fPC1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mPC1", ylab = "fPC1")
ggscatter(Butt, x = "mPC2", y = "fPC2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mPC2", ylab = "fPC2")
malePCs <- ggplot(data = All, aes(x = mPC1, y = mPC2, label = rownames(All))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Male PC1vPC2") + geom_point(aes(colour=All$Strata, size = 4)) 
malePCs + theme_classic() 


M_M <-gls(mPC1~Migration, data=Butt, correlation=Pagel, method = "ML")#Almost, .0523
M_Sa <-gls(mPC1~StrataA, data=Butt, correlation=Pagel) #Under*
M_Sb <-gls(mPC1~Strata, data=Butt, correlation=Pagel)#NOT
M_H <-gls(mPC1~Habitat, data=Butt, correlation=Pagel, method = "ML")
M_F <-gls(mPC1~Forest.Dependency, data=Butt, correlation=Pagel, method = "ML") #NOT
fuck <- M_F <-gls(mPC1~FUCK, data=Butt, correlation=Pagel, method = "ML") 
M_SFa <-gls(mPC1~StrataA+Forest.Dependency, data=Butt, correlation=Pagel) #Under
M_SFb <-gls(mPC1~Strata+Forest.Dependency, data=Butt, correlation=Pagel) #Not
M_SHa <-gls(mPC1~StrataA+Habitat, data=Butt, correlation=Pagel) # Open habitat very, Terriestrial slightly 
M_SHb <-gls(mPC1~Strata+Habitat, data=Butt, correlation=Pagel)#Open habitat very, Open strata slightly
M_MSa <-gls(mPC1~StrataA+Migration, data=Butt, correlation=Pagel)#Understory
M_MSb <-gls(mPC1~Strata+Migration, data=Butt, correlation=Pagel)#NOT
M_HM <-gls(mPC1~Habitat+Migration, data=Butt, correlation=Pagel, method = "ML")#NOT
M_FM <-gls(mPC1~Forest.Dependency+Migration, data=Butt, correlation=Pagel, method = "ML")#NOT
## PC2
M_M <-gls(mPC2~Migration, data=Butt, correlation=Pagel,method = "ML")#NOT
M_Sa <-gls(mPC2~StrataA, data=Butt, correlation=Pagel)#NOT
M_Sb <-gls(mPC2~Strata, data=Butt, correlation=Pagel)#NOT
M_H <-gls(mPC2~Habitat, data=Butt, correlation=Pagel, method = "ML")
M_F <-gls(mPC2~Forest.Dependency, data=Butt, correlation=Pagel, method = "ML") # mEDIUM FOREST SIG
M_SFa <-gls(mPC2~StrataA+Forest.Dependency, data=Butt, correlation=Pagel)# mEDIUM FOREST SIG
M_SFb <-gls(mPC2~Strata+Forest.Dependency, data=Butt, correlation=Pagel)#NOT
M_MSa <-gls(mPC2~StrataA+Migration, data=Butt, correlation=Pagel)#not
M_MSb <-gls(mPC2~Strata+Migration, data=Butt, correlation=Pagel)#NOT
M_SHa <-gls(mPC2~StrataA+Habitat, data=Butt, correlation=Pagel) # Understory 
M_SHb <-gls(mPC2~Strata+Habitat, data=Butt, correlation=Pagel,method = "ML")# fit signular 
M_HM <-gls(mPC2~Habitat+Migration, data=Butt, correlation=Pagel, method = "ML")#Open habitat hIGHLY SIGNIFICANT 
M_FM <-gls(mPC2~Forest.Dependency+Migration, data=Butt, correlation=Pagel, method = "ML")#Low and medium hIGHLY SIGNIFICANT 
a<-AIC(M_M,M_S,M_H,M_F,M_SF,M_MS,M_FM)
summary(M_M)
summary(M_Sa)
summary(M_Sb)
summary(M_H)## 
summary(M_F)##
summary(M_SFa)
summary(M_SFb)
summary(M_HM)# 
summary(M_MSa)
summary(M_MSb)
summary(M_FM)# 
summary(M_SHa)
summary(M_SHb)

F_M <-gls(fPC1~Migration, data=Butt, correlation=Pagel, method = "ML")#almost
F_Sa <-gls(fPC1~StrataA, data=Butt, correlation=Pagel)#not
F_Sb <-gls(fPC1~Strata, data=Butt, correlation=Pagel)#not
F_H <-gls(fPC1~Habitat, data=Butt, correlation=Pagel, method = "ML")#not
F_F <-gls(fPC1~Forest.Dependency, data=Butt, correlation=Pagel, method = "ML") #not
F_HM <-gls(fPC1~Habitat+Migration, data=Butt, correlation=Pagel, method = "ML")#not
F_FM <-gls(fPC1~Forest.Dependency+Migration, data=Butt, correlation=Pagel, method = "ML")#not
F_SFa <-gls(fPC1~StrataA+Forest.Dependency, data=Butt, correlation=Pagel)#not
F_SFb <-gls(fPC1~Strata+Forest.Dependency, data=Butt, correlation=Pagel)#not
F_MSa <-gls(fPC1~StrataA+Migration, data=Butt, correlation=Pagel)#mIGRATION SIG, UNDER SIG
F_MSb <-gls(fPC1~Strata+Migration, data=Butt, correlation=Pagel)#not
F_SHa <-gls(fPC1~StrataA+Habitat, data=Butt, correlation=Pagel) # Not
F_SHb <-gls(fPC1~Strata+Habitat, data=Butt, correlation=Pagel)# fit singular 
## PC2
F_M <-gls(fPC2~Migration, data=Butt, correlation=Pagel, method = "ML")#not
F_Sa <-gls(fPC2~StrataA, data=Butt, correlation=Pagel)#not
F_Sb <-gls(fPC2~Strata, data=Butt, correlation=Pagel)#Sig
F_H <-gls(fPC2~Habitat, data=Butt, correlation=Pagel, method = "ML")#not
F_F <-gls(fPC2~Forest.Dependency, data=Butt, correlation=Pagel, method = "ML") # Best
F_HM <-gls(fPC2~Habitat+Migration, data=Butt, correlation=Pagel, method = "ML")#not
F_FM <-gls(fPC2~Forest.Dependency+Migration, data=Butt, correlation=Pagel, method = "ML")#not
F_SFa <-gls(fPC2~StrataA+Forest.Dependency, data=Butt, correlation=Pagel)#not
F_SFb <-gls(fPC2~Strata+Forest.Dependency, data=Butt, correlation=Pagel)#not
F_MSa <-gls(fPC2~StrataA+Migration, data=Butt, correlation=Pagel)#mIGRATION SIG, UNDER SIG
F_MSb <-gls(fPC2~Strata+Migration, data=Butt, correlation=Pagel)#mIGRATION SIG, UNDER SIG, u/c SIG
F_SHa <-gls(fPC2~StrataA+Habitat, data=Butt, correlation=Pagel) #
F_SHb <-gls(fPC2~Strata+Habitat, data=Butt, correlation=Pagel)#

summary(F_M)
summary(F_Sa)
summary(F_Sb)
summary(F_H)## 
summary(F_F)##
summary(F_HS)#
summary(F_HM)# 
summary(F_MH)# 
summary(F_SFa)
summary(F_SFb)
summary(F_MSa)
summary(F_MSb)
summary(F_SHa)
summary(F_SHb)
#################################
capture.output(summary( a1), file = "p-values/a1.txt")
capture.output(summary( a2), file = "p-values/a2.txt")
capture.output(summary( a3), file = "p-values/a3.txt")
capture.output(summary( a4a), file = "p-values/a4a.txt")
capture.output(summary( a5a), file = "p-values/a5a.txt")
capture.output(summary( a6a), file = "p-values/a6a.txt")
capture.output(summary( a7a), file = "p-values/a7a.txt")
capture.output(summary( a8), file = "p-values/a8.txt")
capture.output(summary( a9), file = "p-values/a9.txt")
capture.output(summary( a10), file = "p-values/a10.txt")
capture.output(summary( a4), file = "p-values/a4.txt")
capture.output(summary( a5), file = "p-values/a5.txt")
capture.output(summary( a6), file = "p-values/a6.txt")
capture.output(summary( b1), file = "p-values/b1.txt")
capture.output(summary( b2), file = "p-values/b2.txt")
capture.output(summary( b3), file = "p-values/b3.txt")
capture.output(summary( b4a), file = "p-values/b4a.txt")
capture.output(summary( b5a), file = "p-values/b5a.txt")
capture.output(summary( b6a), file = "p-values/b6a.txt")
capture.output(summary( b7a), file = "p-values/b7a.txt")
capture.output(summary( b8), file = "p-values/b8.txt")
capture.output(summary( b9), file = "p-values/b9.txt")
capture.output(summary( b10), file = "p-values/b10.txt")
capture.output(summary( b4), file = "p-values/b4.txt")
capture.output(summary( b5), file = "p-values/b5.txt")
capture.output(summary( b6), file = "p-values/b6.txt")
capture.output(summary( c1), file = "p-values/c1.txt")
capture.output(summary( c2), file = "p-values/c2.txt")
capture.output(summary( c3), file = "p-values/c3.txt")
capture.output(summary( c4a), file = "p-values/c4a.txt")
capture.output(summary( c5a), file = "p-values/c5a.txt")
capture.output(summary( c6a), file = "p-values/c6a.txt")
capture.output(summary( c7a), file = "p-values/c7a.txt")
capture.output(summary( c8), file = "p-values/c8.txt")
capture.output(summary( c9), file = "p-values/c9.txt")
capture.output(summary( c10), file = "p-values/c10.txt")
capture.output(summary( c4), file = "p-values/c4.txt")
capture.output(summary( c5), file = "p-values/c5.txt")
capture.output(summary( c6), file = "p-values/c6.txt")
capture.output(summary( d1), file = "p-values/d1.txt")
capture.output(summary( d2), file = "p-values/d2.txt")
capture.output(summary( d3), file = "p-values/d3.txt")
capture.output(summary( d4a), file = "p-values/d4a.txt")
capture.output(summary( d5a), file = "p-values/d5a.txt")
capture.output(summary( d6a), file = "p-values/d6a.txt")
capture.output(summary( d7a), file = "p-values/d7a.txt")
capture.output(summary( d8), file = "p-values/d8.txt")
capture.output(summary( d9), file = "p-values/d9.txt")
capture.output(summary( d10), file = "p-values/d10.txt")
capture.output(summary( d4), file = "p-values/d4.txt")
capture.output(summary( d5), file = "p-values/d5.txt")
capture.output(summary( d6), file = "p-values/d6.txt")

capture.output(summary( w), file = "p-values/w.txt")
capture.output(summary( x), file = "p-values/x.txt")
capture.output(summary( y), file = "p-values/y.txt")
capture.output(summary( z), file = "p-values/z.txt")
capture.output(summary( w1), file = "p-values/w1.txt")
capture.output(summary( x1), file = "p-values/x1.txt")
capture.output(summary( y1), file = "p-values/y1.txt")
capture.output(summary( z1), file = "p-values/z1.txt")
capture.output(summary( w2), file = "p-values/w2.txt")
capture.output(summary( x2), file = "p-values/x2.txt")
capture.output(summary( y2), file = "p-values/y2.txt")
capture.output(summary( z2), file = "p-values/z2.txt")
