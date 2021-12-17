################## Raw Traits Script   #################
################## How raw traits influence plumage evolution
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

setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis")
#FullCardinal.tree <- read.nexus("Cardinalidae.nexus")
#Cardinal.tree <- read.nexus("Cardinalidae.nexus")
tip <- c("Cyanocompsa_cyanoides","Caryothraustes_canadensis")
tip <- c("Cyanocompsa_cyanoides")

Card.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- drop.tip(Card.tree,tip)

This <- read.csv("New.csv")
rownames(All) <- All[,1,]
#All<- All[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15)]
All<- All[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)]
All<-All[Card.tree$tip.label, ] 


############################## Raw Traits
WPTCS <- read.csv("WPTCS.csv")
mWPTCS <- subset(WPTCS, Sex == "Male")
fWPTCS <- subset(WPTCS, Sex == "Female")
#Make species rownames and remove extra cloumn
rownames(mWPTCS) <- mWPTCS[,1,]
mWPTCS<- mWPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
mWPTCS<-mWPTCS[Card.tree$tip.label, ] 
rownames(fWPTCS) <- fWPTCS[,1,]
fWPTCS<- fWPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
fWPTCS<-fWPTCS[Card.tree$tip.label, ] 

#Log-transform 
log.f_WPTSC <- log(fWPTCS[,1:10])
log.m_WPTSC <- log(mWPTCS[,1:10])

#### male PCA 
MPCA <- phyl.pca(Card.tree, mWPTCS, method = "lambda", mode = "cor")
summary(MPCA)#Importance of Components
print(MPCA) #Standard Deviations
m_score <- as.data.frame(MPCA$S) ## Get scores
m_loadings <- as.data.frame(MPCA$L) # Loadings
m_loadings <- m_loadings*-1
#### female PCA 
FPCA <- phyl.pca(Card.tree, fWPTCS, method = "lambda", mode = "cor")
summary(FPCA)#Importance of Components
print(FPCA) #Standard Deviations
f_score <- as.data.frame(FPCA$S)## Get scores
f_loadings <- as.data.frame(FPCA$L) # Loadings
f_loadings <- f_loadings*-1
male <- ggplot(data = This, aes(x = mPC1, y = mPC2, label = rownames(This))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Male PC1vPC2") + geom_point(aes (size = 4)) #(colour=All$Genus,size = 4)) 
male + theme_classic() 

female <- ggplot(data = This, aes(x = fPC1A, y =fPC2A, label = rownames(This))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Female PC1vPC2") + geom_point(aes (size = 4)) #(colour=All$Genus,size = 4)) 
female + theme_classic() 

male <- ggplot(data = This, aes(x = mPC1, y = mPC2, label = rownames(This))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Male PC1vPC2") + geom_point(aes (colour=This$StrataA,size = 4)) #(colour=All$Genus,size = 4)) 
male + theme_classic() 

female <- ggplot(data = This, aes(x = fPC1A, y = fPC2A, label = rownames(This))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Female PC1vPC2") + geom_point(aes (colour=This$Migration,size = 4)) #(colour=All$Genus,size = 4)) 
female + theme_classic() 

male <- ggplot(data = m_score, aes(x = PC1, y = PC2, label = rownames(m_score))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Female PC1vPC2") + geom_point(aes (size = 4)) #(colour=All$Genus,size = 4)) 
male + theme_classic() 

### PGLS carper ##################
####### Caper 

FullCardinal.tree <- read.nexus("Cardinalidae.nexus")
Cardinal.tree <- read.nexus("Cardinalidae.nexus")
#tip <- c("Cyanocompsa_cyanoides","Amaurospiza_concolor")
#Cardinal.tree <- drop.tip(FullCardinal.tree,tip)
Butt <- read.csv("ALL.csv")
WPTCS <- read.csv("WPTCS.csv")
mWPTCS <- subset(WPTCS, Sex == "Male")
fWPTCS <- subset(WPTCS, Sex == "Female")

rownames(mWPTCS) <- mWPTCS[,1,]
mWPTCS<- mWPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
mWPTCS<-mWPTCS[Cardinal.tree$tip.label, ] 
rownames(fWPTCS) <- fWPTCS[,1,]
fWPTCS<- fWPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
fWPTCS<-fWPTCS[Cardinal.tree$tip.label, ] 
log.f_WPTSC <- log(fWPTCS[,1:10])
log.m_WPTSC <- log(mWPTCS[,1:10])
#####
fBrill <- log.f_WPTSC[,"AvgBrill"]
fChroma <- log.f_WPTSC[,"AvgChroma"]
fASpan <- log.f_WPTSC[,"AvgSpan"]
fMSpan <- log.f_WPTSC[,"MaxSpan"]
fVol <- log.f_WPTSC[,"Volume"]
fAHue <- log.f_WPTSC[,"AvgHueDisp"]
fMHue <- log.f_WPTSC[,"MaxHueDisp"]

mBrill <- log.m_WPTSC[,"AvgBrill"]
mChroma <- log.m_WPTSC[,"AvgChroma"]
mASpan <- log.m_WPTSC[,"AvgSpan"]
mMSpan <- log.m_WPTSC[,"MaxSpan"]
mVol <- log.m_WPTSC[,"Volume"]
mAHue <- log.m_WPTSC[,"AvgHueDisp"]
mMHue <- log.m_WPTSC[,"MaxHueDisp"]

names(fBrill) <- names(fChroma) <- names(fASpan)<- names(fMSpan) <- names(fVol) <- names(fAHue) <- names(fMHue) <- rownames(log.f_WPTSC)
names(mBrill) <- names(mChroma) <- names(fASpan)<- names(fMSpan)<-names(mVol)<- names(mAHue) <- names(mMHue) <- rownames(log.m_WPTSC)
Butt <- cbind(Butt,fBrill,mBrill, fChroma, mChroma, mASpan, fASpan, mMSpan, fMSpan, mVol, fVol,mAHue, mMHue, fAHue, fMHue)


comp.data<-comparative.data(Cardinal.tree, Butt, names.col = "ï..Species", vcv.dim=2, warn.dropped=TRUE)

PC1<-pgls(mPC1~fPC1, data=comp.data,lambda =  "ML")
PC2<-pgls(mPC2~fPC2, data=comp.data,lambda =  "ML")
Brill<-pgls(mBrill~fBrill, data=comp.data,lambda =  "ML") # all plots looks good 
Chroma<-pgls(mChroma~fChroma, data=comp.data,lambda =  "ML")# all plots looks good 
Volume<-pgls(mVol~fVol, data=comp.data,lambda =  "ML")
AvgSpan<-pgls(mASpan~fASpan, data=comp.data,lambda =  "ML")
MaxSpan<-pgls(mMSpan~fMSpan, data=comp.data,lambda =  "ML") # Use maxSpan over average span
AvgHue<-pgls(mAHue~fAHue, data=comp.data,lambda =  "ML") ### All plots look good
MaxHue<-pgls(mMHue~fMHue, data=comp.data,lambda =  "ML") ## Could be better 
plot(PC2)
summary(PC1)




#### Brillance

#Migration
m1<-pgls(mBrill~Migration, data=comp.data,lambda = "ML") ### Not
f1<-pgls(fBrill~Migration, data=comp.data,lambda = "ML")## Not
summary(m1)
summary(f1)


#Forest
m2<-pgls(mBrill~Forest.Dependency, data=comp.data,lambda = "ML") ### Not
f2<-pgls(fBrill~Forest.Dependency, data=comp.data,lambda = "ML")## Medium


#Habitat
m3<-pgls(mBrill~Habitat, data=comp.data,lambda = "ML") ### Not  ### @nd best model selection 
f3<-pgls(fBrill~Habitat, data=comp.data,lambda = "ML")## Not
summary(m3)
summary(f3)



## Strata  
m4a<-pgls(mBrill~StrataA, data=comp.data,lambda = "ML")## Under** sig ## 3rd best model selection 
m4<-pgls(mBrill~Strata, data=comp.data,lambda = "ML")## Under* sig
f4a<-pgls(fBrill~StrataA, data=comp.data,lambda = "ML")## None
f4<-pgls(fBrill~Strata, data=comp.data,lambda = "ML")## None


## Strata plus Migration
m5a<-pgls(mBrill~StrataA+Migration, data=comp.data,lambda = "ML")## Under sig**
m5<-pgls(mBrill~Strata+Migration, data=comp.data,lambda = "ML")## None
f5a<-pgls(fBrill~StrataA+Migration, data=comp.data,lambda = "ML")## None
f5<-pgls(fBrill~Strata+Migration, data=comp.data,lambda = "ML")## None

## Strata plus Forest
m6a<-pgls(mBrill~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Under** highly sig ### Best model selection
m6<-pgls(mBrill~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Under highly sig
f6a<-pgls(fBrill~StrataA+Forest.Dependency, data=comp.data,lambda = "ML")## Medium sig
f6<-pgls(fBrill~Strata+Forest.Dependency, data=comp.data,lambda = "ML")## Medium**, no forest*

summary(m6a)
## Strata plus Habitat ==> StrataA good, Strata is singular  
m7a<-pgls(mBrill~StrataA+Habitat, data=comp.data,kappa = "ML")## terrestrial*, U/C, open habitat**
m7<-pgls(mBrill~Strata+Habitat, data=comp.data,delta = "ML")## Singular
f7a<-pgls(fBrill~StrataA+Habitat, data=comp.data,lambda = "ML")## None
f7<-pgls(fBrill~Strata+Habitat, data=comp.data,lambda = "ML")## Singular

# Habitat + Migration 
m8<-pgls(mBrill~Habitat+Migration, data=comp.data,lambda = "ML") ### no
f8<-pgls(fBrill~Habitat+Migration, data=comp.data,lambda = "ML")## no  #### Best model selection 

# Habitat + Forest 
m9<-pgls(mBrill~Forest.Dependency+Habitat, data=comp.data,lambda = "ML") ### no
f9<-pgls(fBrill~Forest.Dependency+Habitat, data=comp.data,lambda = "ML")## Medium **  ### 2nd best model selection 

# Forest + Migration 
m10<-pgls(mBrill~Forest.Dependency+Migration, data=comp.data,lambda = "ML") ### no
f10<-pgls(fBrill~Forest.Dependency+Migration, data=comp.data,lambda = "ML")## Medium *

#Forest vs nonforest
m11<-pgls(mBrill~FUCK, data=comp.data,lambda = "ML") ### Not
f11<-pgls(fBrill~FUCK, data=comp.data,lambda = "ML")## Medium

#Forest vs nonforest
m12<-pgls(mBrill~FUCK+StrataA, data=comp.data,lambda = "ML") ### Not
f12<-pgls(fBrill~FUCK+StrataA, data=comp.data,lambda = "ML")## Medium
summary(m12)

m13<-pgls(mBrill~FUCK+Migration, data=comp.data,lambda = "ML") ### Not
f13<-pgls(fBrill~FUCK+Migration, data=comp.data,lambda = "ML")## Medium

m14<-pgls(mBrill~FUCK+Forest.Dependency, data=comp.data,lambda = "ML") ### Not
f14<-pgls(fBrill~FUCK+Forest.Dependency, data=comp.data,lambda = "ML")## Medium

m15<-pgls(mBrill~FUCK+Migration+StrataA, data=comp.data,lambda = "ML") ### Not
f15<-pgls(fBrill~FUCK+Migration+StrataA, data=comp.data,lambda = "ML")## Medium

m16 <-pgls(mBrill~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
m17 <-pgls(mBrill~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 
f16 <-pgls(fBrill~StrataA+Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
f17 <-pgls(fBrill~StrataA+Habitat+Migration, data=comp.data,lambda = "ML")# 


### AICc model selection 
#1/a1)Migration, 2/a2)Forest.Dep, 3/a3) Habitat, 4/a4a)StrataA, 5/a5a)StrataA+Migration, 6/a6a)StrataA+Forest.Dep, 
# 7/a8)Habitat+Migration  8/a9) Forest.Dep+Habitat  9/a10)Forest.Dep+Migration
# 10/a4) Strata   11/a5) Strata+Migration   12/a6)Strata+Forest.Dep  
# 13/w) nonforest 14/w1) nonforest+StrataA 15/w2)nonforest+Migration  
#16/w4)nonforest+Migration+StrataA, 17/mPC1_SFM) StrataA+Forest.Dep+Habitat 
mBrillaicc <- cbind(m1$aicc,m2$aicc,m3$aicc,m4a$aicc,m5a$aicc,m6a$aicc,m8$aicc,m9$aicc,m10$aicc,m4$aicc,
                    m5$aicc,m6$aicc,m11$aicc,m12$aicc,m13$aicc,m15$aicc,m16$aicc,m17$aicc) 

fBrillaicc <- cbind(f1$aicc,f2$aicc,f3$aicc,f4a$aicc,f5a$aicc,f6a$aicc,f8$aicc,f9$aicc,f10$aicc,f4$aicc,
                    f5$aicc,f6$aicc,f11$aicc,f12$aicc,f13$aicc,f15$aicc,f16$aicc,f17$aicc) 
mBrillaicc
fBrillaicc

Female <- ggplot(data = All, aes(x = fPC1, y = fPC2, label = rownames(All))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") + 
  ggtitle("Female PC1vPC2") + geom_point(aes (colour=All$Migration,size = 4)) #(colour=All$Genus,size = 4)) 
Female + theme_classic() 

####################################
## Save results as trext file
capture.output(summary( fBrill_M), file = "p-values/Brillance/fBrill_M.txt")
capture.output(summary( fBrill_Sa), file = "p-values/Brillance/fBrill_Sa.txt")
capture.output(summary( fBrill_Sb), file = "p-values/Brillance/fBrill_Sb.txt")
capture.output(summary( fBrill_H), file = "p-values/Brillance/fBrill_H.txt")
capture.output(summary( fBrill_F), file = "p-values/Brillance/fBrill_F.txt")
capture.output(summary( fBrill_SFa), file = "p-values/Brillance/fBrill_SFa.txt")
capture.output(summary( fBrill_SFb), file = "p-values/Brillance/fBrill_SFb.txt")
capture.output(summary( fBrill_SHa), file = "p-values/Brillance/fBrill_SHa.txt")


####### Making Figures ####################

fBrill <- log.f_WPTSC[,"AvgBrill"]
fChroma <- log.f_WPTSC[,"AvgChroma"]
fASpan <- log.f_WPTSC[,"AvgSpan"]
fMSpan <- log.f_WPTSC[,"MaxSpan"]
fVol <- log.f_WPTSC[,"Volume"]
mBrill <- log.m_WPTSC[,"AvgBrill"]
mChroma <- log.m_WPTSC[,"AvgChroma"]
mASpan <- log.m_WPTSC[,"AvgSpan"]
mMSpan <- log.m_WPTSC[,"MaxSpan"]
mVol <- log.m_WPTSC[,"Volume"]
mAHue <- log.m_WPTSC[,"AvgHueDisp"]
mMHue <- log.m_WPTSC[,"MaxHueDisp"]
fAHue <- log.f_WPTSC[,"AvgHueDisp"]
fMHue <- log.f_WPTSC[,"MaxHueDisp"]
names(fBrill) <- names(fChroma) <- names(fASpan) <- names(fMSpan) <- names(fVol) <- names(fAHue) <- names(fMHue) <- rownames(log.f_WPTSC)
names(mBrill) <- names(mChroma) <- names(mASpan)<- names(mMSpan) <-names(mVol)<- names(mAHue) <- names(mMHue) <- rownames(log.m_WPTSC)
All <- cbind(All,fBrill,mBrill, fChroma, mChroma, mASpan, fASpan,mMSpan, fMSpan, mVol, fVol,mAHue, mMHue, fAHue, fMHue)

bm<-corBrownian(1, Card.tree) # Random walk 
Pagel<-corPagel(1, Card.tree) # Allows BM to devate slightly
OU<-corMartins(1, Card.tree) # Evolving under stabilizing selection, or rubber brand
EB <- corBlomberg(1, Card.tree) # (EB) Assume BM but rates accelerate (g<1) or decelerate (g>1) through time. Evolve Under adaptive Radiation? 

# PC1
M_F <-gls(mPC1~fPC1, data=All, correlation=Pagel)
summary(M_F) # Correlated ***
plot(mPC1~fPC1, data = All)
abline(M_F)
M_F <-gls(mPC2~fPC2, data=All, correlation=Pagel)
summary(M_F) # Correlated ***
plot(mPC2~fPC2, data = All)
abline(M_F)
# Brillance
M1_F1 <-gls(mBrill~fBrill, data=All, correlation=Pagel)
summary(M1_F1) # Correlated ***
plot(mBrill~fBrill, data = All)
abline(M1_F1)
###### Chroma
M2_F2 <-gls(mChroma~fChroma, data=All, correlation=Pagel)
summary(M2_F2) # Correlated ***
plot(mChroma~fChroma, data = All)
abline(M2_F2)
# Vol
M3_F3 <-gls(mVol~fVol, data=All, correlation=Pagel)
summary(M3_F3) # Correlated **
plot(mVol~fVol, data = All)
abline(M3_F3)
#AvgSpan
M4_F4 <-gls(mASpan~fASpan, data=All, correlation=Pagel)
summary(M4_F4) # Correlated **
plot(mASpan~fASpan, data = All)
abline(M4_F4)
#MaxSpan
M5_F5 <-gls(mMSpan~fMSpan, data=All, correlation=Pagel)
summary(M5_F5) # Correlated *
plot(mMSpan~fMSpan, data = All)
abline(M5_F5)
#Average Hue
M6_F6 <-gls(mAHue~fAHue, data=All, correlation=Pagel)
summary(M6_F6) # Not
plot(mAHue~fAHue, data = All)
abline(M6_F6)
#Max Hue
M7_F7 <-gls(mMHue~fMHue, data=All, correlation=Pagel)
summary(M7_F7) # ***
plot(mMhue~fMhue, data = All)
abline(M7_F7)

ggscatter(All, x = "mBrill", y = "fBrill", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mBrill", ylab = "fBrill")
ggscatter(All, x = "mChroma", y = "fChroma", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mChroma", ylab = "fChroma")
ggscatter(All, x = "mSpan", y = "fSpan", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mSpan", ylab = "fSpan")
ggscatter(All, x = "mVol", y = "fVol", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mVol", ylab = "fVol")
ggscatter(All, x = "mAHue", y = "fAHue", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mAhue", ylab = "fAhue")


####### Other traits #########################
#Hue
mAhue_M <-pgls(mAHue~Migration, data=comp.data,lambda = "ML") 
mAhue_Sa <-pgls(mAHue~StrataA, data=comp.data,lambda = "ML") # Understory ***
mAhue_Sb <-pgls(mAHue~Strata, data=comp.data,lambda = "ML")#Understory ***
mAhue_H <-pgls(mAHue~Habitat, data=comp.data,lambda = "ML") # Not
mAhue_F <-pgls(mAHue~Forest.Dependency, data=comp.data,lambda = "ML") ## Not
mAhue_SFa <-pgls(mAHue~StrataA+Forest.Dependency, data=comp.data,lambda = "ML") ## Understory 
mAhue_SFb <-pgls(mAHue~Strata+Forest.Dependency, data=comp.data,lambda = "ML") ### Understory 
mAhue_SHa <-pgls(mAHue~StrataA+Habitat, data=comp.data,lambda = "ML") # ### Understory 
mAhue_SHb <-pgls(mAHue~Strata+Habitat, data=comp.data,lambda = "ML")#### Understory 
mAhue_MSa <-pgls(mAHue~StrataA+Migration, data=comp.data,lambda = "ML")#### Understory 
mAhue_MSb <-pgls(mAHue~Strata+Migration, data=comp.data,lambda = "ML")#### Understory 
mAhue_HM <-pgls(mAHue~Habitat+Migration, data=comp.data,lambda = "ML")# Migration yes
mAhue_FM <-pgls(mAHue~Forest.Dependency+Migration, data=comp.data,lambda = "ML")# Migration
summary(mAhue_M)# Not
summary(mAhue_Sa)# Not
summary(mAhue_Sb)# Not
summary(mAhue_H)## # Not
summary(mAhue_F)### # Medium Forest*
summary(mAhue_SFa)### Medium Forest*
summary(mAhue_SFb)### Medium Forest*
summary(mAhue_HM)# # Not
summary(mAhue_MSa)# Not
summary(mAhue_MSb)# Not
summary(mAhue_FM)# ### Medium Forest*
summary(mAhue_SHa)# Not
summary(mAhue_SHb)# Not
# 3) Female AHue
fAhue_M <-pgls(fAHue~Migration, data=comp.data,lambda = "ML")# 
fAhue_Sa <-pgls(fAHue~StrataA, data=comp.data,lambda = "ML") # 
fAhue_Sb <-pgls(fAHue~Strata, data=comp.data,lambda = "ML")#
fAhue_H <-pgls(fAHue~Habitat, data=comp.data,lambda = "ML") # 
fAhue_F <-pgls(fAHue~Forest.Dependency, data=comp.data,lambda = "ML") ## 
fAhue_SFa <-pgls(fAHue~StrataA+Forest.Dependency, data=comp.data,lambda = "ML") ## 
fAhue_SFb <-pgls(fAHue~Strata+Forest.Dependency, data=comp.data,lambda = "ML") ### 
fAhue_SHa <-pgls(fAHue~StrataA+Habitat, data=comp.data,lambda = "ML") # ### 
fAhue_SHb <-pgls(fAHue~Strata+Habitat, data=comp.data,lambda = "ML")#### 
fAhue_MSa <-pgls(fAHue~StrataA+Migration, data=comp.data,lambda = "ML")#### 
fAhue_MSb <-pgls(fAHue~Strata+Migration, data=comp.data,lambda = "ML")#### 
fAhue_HM <-pgls(fAHue~Habitat+Migration, data=comp.data,lambda = "ML")# 
fAhue_FM <-pgls(fAHue~Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
summary(fAhue_M)## Not
summary(fAhue_Sa)## Not
summary(fAhue_Sb)## Not
summary(fAhue_H)## Not
summary(fAhue_F)#### Not
summary(fAhue_HS)### Not
summary(fAhue_MH)# ## Not
summary(fAhue_SFa)## Not
summary(fAhue_SFb)## Not
summary(fAhue_SHa)## Not
summary(fAhue_SHb)## Not
summary(fAhue_MSa)### Understory, Sedintary
summary(fAhue_MSb)## Not
summary(fAhue_HM)### Not
summary(fAhue_FM)# ## Not


#Male Volume
mVol_M <-pgls(mVol~Migration, data=comp.data,lambda = "ML") 
mVol_Sa <-pgls(mVol~StrataA, data=comp.data,lambda = "ML") # 
mVol_Sb <-pgls(mVol~Strata, data=comp.data,lambda = "ML")#
mVol_H <-pgls(mVol~Habitat, data=comp.data,lambda = "ML") # Not
mVol_F <-pgls(mVol~Forest.Dependency, data=comp.data,lambda = "ML") ## Not
mVol_SFa <-pgls(mVol~StrataA+Forest.Dependency, data=comp.data,lambda = "ML") ##  
mVol_SFb <-pgls(mVol~Strata+Forest.Dependency, data=comp.data,lambda = "ML") ### 
mVol_SHa <-pgls(mVol~StrataA+Habitat, data=comp.data,lambda = "ML") #
mVol_SHb <-pgls(mVol~Strata+Habitat, data=comp.data,lambda = "ML")#### 
mVol_MSa <-pgls(mVol~StrataA+Migration, data=comp.data,lambda = "ML")#### 
mVol_MSb <-pgls(mVol~Strata+Migration, data=comp.data,lambda = "ML")#### 
mVol_HM <-pgls(mVol~Habitat+Migration, data=comp.data,lambda = "ML")# 
mVol_FM <-pgls(mVol~Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
summary(mVol_M)# Almost
summary(mVol_Sa)# Not
summary(mVol_Sb)# Not
summary(mVol_H)## # Not
summary(mVol_F)### # 
summary(mVol_SFa)### 
summary(mVol_SFb)### 
summary(mVol_HM)# # Migration
summary(mVol_MSa)# Not
summary(mVol_MSb)# Not
summary(mVol_FM)# ### Not
summary(mVol_SHa)# OPen Habitat Highly correlated with Volume 
summary(mVol_SHb)#  
# 3) Female Volume
fVol_M <-pgls(fVol~Migration, data=comp.data,lambda = "ML")# 
fVol_Sa <-pgls(fVol~StrataA, data=comp.data,lambda = "ML") # 
fVol_Sb <-pgls(fVol~Strata, data=comp.data,lambda = "ML")#
fVol_H <-pgls(fVol~Habitat, data=comp.data,lambda = "ML") # 
fVol_F <-pgls(fVol~Forest.Dependency, data=comp.data,lambda = "ML") ## 
fVol_SFa <-pgls(fVol~StrataA+Forest.Dependency, data=comp.data,lambda = "ML") ## 
fVol_SFb <-pgls(fVol~Strata+Forest.Dependency, data=comp.data,lambda = "ML") ### 
fVol_SHa <-pgls(fVol~StrataA+Habitat, data=comp.data,lambda = "ML") # ### 
fVol_SHb <-pgls(fVol~Strata+Habitat, data=comp.data,lambda = "ML")#### 
fVol_MSa <-pgls(fVol~StrataA+Migration, data=comp.data,lambda = "ML")#### 
fVol_MSb <-pgls(fVol~Strata+Migration, data=comp.data,lambda = "ML")#### 
fVol_HM <-pgls(fVol~Habitat+Migration, data=comp.data,lambda = "ML")# 
fVol_FM <-pgls(fVol~Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
summary(fVol_M)## Not
summary(fVol_Sa)## Understoryy*
summary(fVol_Sb)## Not
summary(fVol_H)## Not
summary(fVol_F)#### Not
summary(fVol_SFa)## understory*
summary(fVol_SFb)## Not
summary(fVol_SHa)## understory*
summary(fVol_SHb)## Not
summary(fVol_MSa)### Understory*** Sedintary***
summary(fVol_MSb)## Understory* Sedintary*
summary(fVol_HM)### Not
summary(fVol_FM)# ## Not

###### Chroma

#Male Volume
mChroma_M <-pgls(mChroma~Migration, data=comp.data,lambda = "ML") 
mChroma_Sa <-pgls(mChroma~StrataA, data=comp.data,lambda = "ML") # 
mChroma_Sb <-pgls(mChroma~Strata, data=comp.data,lambda = "ML")#
mChroma_H <-pgls(mChroma~Habitat, data=comp.data,lambda = "ML") # Not
mChroma_F <-pgls(mChroma~Forest.Dependency, data=comp.data,lambda = "ML") ## Not
mChroma_SFa <-pgls(mChroma~StrataA+Forest.Dependency, data=comp.data,lambda = "ML") ##  
mChroma_SFb <-pgls(mChroma~Strata+Forest.Dependency, data=comp.data,lambda = "ML") ### 
mChroma_SHa <-pgls(mChroma~StrataA+Habitat, data=comp.data,lambda = "ML") #
mChroma_SHb <-pgls(mChroma~Strata+Habitat, data=comp.data,lambda = "ML")#### 
mChroma_MSa <-pgls(mChroma~StrataA+Migration, data=comp.data,lambda = "ML")#### 
mChroma_MSb <-pgls(mChroma~Strata+Migration, data=comp.data,lambda = "ML")#### 
mChroma_HM <-pgls(mChroma~Habitat+Migration, data=comp.data,lambda = "ML")# 
mChroma_FM <-pgls(mChroma~Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
summary(mChroma_M)# Almost
summary(mChroma_Sa)# Not
summary(mChroma_Sb)# Not
summary(mChroma_H)## # Almost open
summary(mChroma_F)### # 
summary(mChroma_SFa)### 
summary(mChroma_SFb)### 
summary(mChroma_HM)# # Migration
summary(mChroma_MSa)# Not
summary(mChroma_MSb)# Not
summary(mChroma_FM)# ### Not
summary(mChroma_SHa)# OPen Habitat Highly correlated with Chroma, Terrestrial  
 
# 3) Female Volume
fChroma_M <-pgls(fChroma~Migration, data=comp.data,lambda = "ML")# 
fChroma_Sa <-pgls(fChroma~StrataA, data=comp.data,lambda = "ML") # 
fChroma_Sb <-pgls(fChroma~Strata, data=comp.data,lambda = "ML")#
fChroma_H <-pgls(fChroma~Habitat, data=comp.data,lambda = "ML") # 
fChroma_F <-pgls(fChroma~Forest.Dependency, data=comp.data,lambda = "ML") ## 
fChroma_SFa <-pgls(fChroma~StrataA+Forest.Dependency, data=comp.data,lambda = "ML") ## 
fChroma_SFb <-pgls(fChroma~Strata+Forest.Dependency, data=comp.data,lambda = "ML") ### 
fChroma_SHa <-pgls(fChroma~StrataA+Habitat, data=comp.data,lambda = "ML") # ### 
fChroma_SHb <-pgls(fChroma~Strata+Habitat, data=comp.data,lambda = "ML")#### 
fChroma_MSa <-pgls(fChroma~StrataA+Migration, data=comp.data,lambda = "ML")#### 
fChroma_MSb <-pgls(fChroma~Strata+Migration, data=comp.data,lambda = "ML")#### 
fChroma_HM <-pgls(fChroma~Habitat+Migration, data=comp.data,lambda = "ML")# 
fChroma_FM <-pgls(fChroma~Forest.Dependency+Migration, data=comp.data,lambda = "ML")# 
summary(fChroma_M)## Not
summary(fChroma_Sa)## not
summary(fChroma_Sb)## Not
summary(fChroma_H)## Not
summary(fChroma_F)#### Not
summary(fChroma_SFa)## understory*
summary(fChroma_SFb)## Not
summary(fChroma_SHa)## understory*
summary(fChroma_SHb)## Not
summary(fChroma_MSa)### Understory*** Sedintary***
summary(fChroma_MSb)## Understory* Sedintary*
summary(fChroma_HM)### Not
summary(fChroma_FM)# ## Not
