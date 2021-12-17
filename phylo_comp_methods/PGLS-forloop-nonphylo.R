#### PGLS for-loop #### 
library(ape)
library(phytools)
library(Rcpp)
library(geiger)
library(RPANDA)
library(caper)
library(MCMCglmm)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(MuMIn)

########################################################################

setwd("C:/Users/bfsco/Desktop/Masters Research/Chapter 2/analysis/WPTSC/MT_tree_analysis")

# Read in Tree and data 
Card.tree <- read.nexus("FullMCC.tree.nexus")
plot(Card.tree)
#newer <- read.csv("IUCN/new.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') 
data <- read.csv("Chap2data.nonphylo.csv", header = TRUE) 

rownames(data) <- data[,1,]
Brill <- log(data[c(13,14)])
data <- data[c(-13,-14)]
df <- cbind(data,Brill)


# need to drop this species from all analysis 
# "Periporphyrus_erythromelas"

#Align species name to tree
df<-df[Card.tree$tip.label, ]
comp.data<-comparative.data(Card.tree, data, names.col = "phylo", vcv=TRUE, warn.dropped=TRUE)

# Comparing male and female plumage 
a <- pgls(data=comp.data,fPC1~mPC1,lambda="ML")
summary(a)
with(data, plot(fPC1~mPC1, xlab = "fPC1", ylab = "mPC1",
                main = "fPC1 vs \ mPC1"))
abline(a)

b <- pgls(data=comp.data,fPC2~mPC2,lambda="ML")
summary(b)
b$aicc

with(data, plot(fPC2~mPC2, xlab = "fPC2", ylab = "mPC2",
                main = "fPC2 vs \ mPC2"))
abline(b)


c <- pgls(fPC1~Habitat+Migration,data=comp.data,lambda="ML")
summary(c)


Habitat

d <- pgls(mAvgBrill~nonforest + Strata+Migration, data=comp.data,lambda="ML")
summary(d)
d$aicc


###################
FvMmods<-list()
FvMmods[[1]]<-pgls(data=comp.data,fPC1~mPC1,lambda="ML")
FvMmods[[2]]<-pgls(data=comp.data,fPC2~mPC2,lambda="ML")
FvMmods[[3]]<-pgls(data=comp.data,fPC3~mPC3,lambda="ML")

FvMmodsAIC<-matrix(nrow=3,ncol=5)
for (i in 1:3){
  FvMmodsAIC[i,2]<-FvMmods[[i]]$aicc
  FvMmodsAIC[i,3]<-summary(FvMmods[[i]])$adj.r.squared
  FvMmodsAIC[i,4]<-pf(summary(FvMmods[[i]])$fstatistic[1],summary(FvMmods[[i]])$fstatistic[2],summary(FvMmods[[i]])$fstatistic[3],lower.tail=F)
  FvMmodsAIC[i,5]<-FvMmods[[i]]$param.CI$lambda$opt
  FvMmodsAIC[i,1]<-paste("FvM_Model",i,sep="_")
}

FvMmodsAIC<-data.frame(FvMmodsAIC)
FvMmodsAIC[,2]<-as.numeric(as.character(FvMmodsAIC[,2]))
FvMmodsAIC[,3]<-as.numeric(as.character(FvMmodsAIC[,3]))
FvMmodsAIC[,5]<-as.numeric(as.character(FvMmodsAIC[,5]))

names(FvMmodsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda")

FvM<-FvMmodsAIC[order(FvMmodsAIC$AICc),]
write.csv(FvM,"Results/PGLS.nonphylo/FvM_models.csv")


############################################
#https://github.com/enicurus/warbler.molt.migration/blob/master/mixed_models.R


# Just categorical variables 
mPC1mods<-list()
mPC1mods[[1]]<-pgls(data=comp.data,mPC1~nonforest,lambda="ML")
mPC1mods[[2]]<-pgls(data=comp.data,mPC1~Strata,lambda="ML")
mPC1mods[[3]]<-pgls(data=comp.data,mPC1~Migration,lambda="ML")
mPC1mods[[4]]<-pgls(data=comp.data,mPC1~Forest.Dependency,lambda="ML")
mPC1mods[[5]]<-pgls(data=comp.data,mPC1~nonforest+Strata,lambda="ML")
mPC1mods[[6]]<-pgls(data=comp.data,mPC1~nonforest+Migration,lambda="ML")
mPC1mods[[7]]<-pgls(data=comp.data,mPC1~Forest.Dependency+Strata,lambda="ML")
mPC1mods[[8]]<-pgls(data=comp.data,mPC1~Forest.Dependency+Migration,lambda="ML")
mPC1mods[[9]]<-pgls(data=comp.data,mPC1~Migration+Strata,lambda="ML")
mPC1mods[[10]]<-pgls(data=comp.data,mPC1~Migration+Strata+nonforest,lambda="ML")
mPC1mods[[11]]<-pgls(data=comp.data,mPC1~Migration+Strata+Forest.Dependency,lambda="ML")
mPC1mods[[12]]<-pgls(data=comp.data,mPC1~Habitat,lambda="ML")
mPC1mods[[13]]<-pgls(data=comp.data,mPC1~Habitat+Strata,lambda="ML")
mPC1mods[[14]]<-pgls(data=comp.data,mPC1~Habitat+Migration,lambda="ML")
mPC1mods[[15]]<-pgls(data=comp.data,mPC1~Migration+Strata+Habitat,lambda="ML")

mPC1modsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  mPC1modsAIC[i,2]<-mPC1mods[[i]]$aicc
  mPC1modsAIC[i,3]<-summary(mPC1mods[[i]])$adj.r.squared
  mPC1modsAIC[i,4]<-pf(summary(mPC1mods[[i]])$fstatistic[1],summary(mPC1mods[[i]])$fstatistic[2],summary(mPC1mods[[i]])$fstatistic[3],lower.tail=F)
  mPC1modsAIC[i,5]<-mPC1mods[[i]]$param.CI$lambda$opt
  mPC1modsAIC[i,1]<-paste("mPC1_Model",i,sep="_")
}

mPC1modsAIC<-data.frame(mPC1modsAIC)
mPC1modsAIC[,2]<-as.numeric(as.character(mPC1modsAIC[,2]))
mPC1modsAIC[,3]<-as.numeric(as.character(mPC1modsAIC[,3]))
mPC1modsAIC[,5]<-as.numeric(as.character(mPC1modsAIC[,5]))

names(mPC1modsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda")

mPC1<-mPC1modsAIC[order(mPC1modsAIC$AICc),]
write.csv(mPC1,"Results/PGLS.nonphylo/mPC1_models.csv")

bonfor.mpc1 <- p.adjust(mPC1$`p-value`, method = p.adjust.methods, n = length(mPC1$`p-value`))


##############################
# Just categorical variables 
mPC2mods<-list()
mPC2mods[[1]]<-pgls(data=comp.data,mPC2~nonforest,lambda="ML")
mPC2mods[[2]]<-pgls(data=comp.data,mPC2~Strata,lambda="ML")
mPC2mods[[3]]<-pgls(data=comp.data,mPC2~Migration,lambda="ML")
mPC2mods[[4]]<-pgls(data=comp.data,mPC2~Forest.Dependency,lambda="ML")
mPC2mods[[5]]<-pgls(data=comp.data,mPC2~nonforest+Strata,lambda="ML")
mPC2mods[[6]]<-pgls(data=comp.data,mPC2~nonforest+Migration,lambda="ML")
mPC2mods[[7]]<-pgls(data=comp.data,mPC2~Forest.Dependency+Strata,lambda="ML")
mPC2mods[[8]]<-pgls(data=comp.data,mPC2~Forest.Dependency+Migration,lambda="ML")
mPC2mods[[9]]<-pgls(data=comp.data,mPC2~Migration+Strata,lambda="ML")
mPC2mods[[10]]<-pgls(data=comp.data,mPC2~Migration+Strata+nonforest,lambda="ML")
mPC2mods[[11]]<-pgls(data=comp.data,mPC2~Migration+Strata+Forest.Dependency,lambda="ML")
mPC2mods[[12]]<-pgls(data=comp.data,mPC2~Habitat,lambda="ML")
mPC2mods[[13]]<-pgls(data=comp.data,mPC2~Habitat+Strata,lambda="ML")
mPC2mods[[14]]<-pgls(data=comp.data,mPC2~Habitat+Migration,lambda="ML")
mPC2mods[[15]]<-pgls(data=comp.data,mPC2~Migration+Strata+Habitat,lambda="ML")

mPC2modsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  mPC2modsAIC[i,2]<-mPC2mods[[i]]$aicc
  mPC2modsAIC[i,3]<-summary(mPC2mods[[i]])$adj.r.squared
  mPC2modsAIC[i,4]<-pf(summary(mPC2mods[[i]])$fstatistic[1],summary(mPC2mods[[i]])$fstatistic[2],summary(mPC2mods[[i]])$fstatistic[3],lower.tail=F)
  mPC2modsAIC[i,5]<-mPC2mods[[i]]$param.CI$lambda$opt
  mPC2modsAIC[i,1]<-paste("mPC2_Model",i,sep="_")
}

mPC2modsAIC<-data.frame(mPC2modsAIC)
mPC2modsAIC[,2]<-as.numeric(as.character(mPC2modsAIC[,2]))
mPC2modsAIC[,3]<-as.numeric(as.character(mPC2modsAIC[,3]))
mPC2modsAIC[,5]<-as.numeric(as.character(mPC2modsAIC[,5]))

names(mPC2modsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda")

mPC2<-mPC2modsAIC[order(mPC2modsAIC$AICc),]
write.csv(mPC2,"Results/PGLS.nonphylo/mPC2_models.csv")

################################
# Just categorical variables 
# Just categorical variables 
fPC1mods<-list()
fPC1mods[[1]]<-pgls(data=comp.data,fPC1~nonforest,lambda="ML")
fPC1mods[[2]]<-pgls(data=comp.data,fPC1~Strata,lambda="ML")
fPC1mods[[3]]<-pgls(data=comp.data,fPC1~Migration,lambda="ML")
fPC1mods[[4]]<-pgls(data=comp.data,fPC1~Forest.Dependency,lambda="ML")
fPC1mods[[5]]<-pgls(data=comp.data,fPC1~nonforest+Strata,lambda="ML")
fPC1mods[[6]]<-pgls(data=comp.data,fPC1~nonforest+Migration,lambda="ML")
fPC1mods[[7]]<-pgls(data=comp.data,fPC1~Forest.Dependency+Strata,lambda="ML")
fPC1mods[[8]]<-pgls(data=comp.data,fPC1~Forest.Dependency+Migration,lambda="ML")
fPC1mods[[9]]<-pgls(data=comp.data,fPC1~Migration+Strata,lambda="ML")
fPC1mods[[10]]<-pgls(data=comp.data,fPC1~Migration+Strata+nonforest,lambda="ML")
fPC1mods[[11]]<-pgls(data=comp.data,fPC1~Migration+Strata+Forest.Dependency,lambda="ML")
fPC1mods[[12]]<-pgls(data=comp.data,fPC1~Habitat,lambda="ML")
fPC1mods[[13]]<-pgls(data=comp.data,fPC1~Habitat+Strata,lambda="ML")
fPC1mods[[14]]<-pgls(data=comp.data,fPC1~Habitat+Migration,lambda="ML")
fPC1mods[[15]]<-pgls(data=comp.data,fPC1~Migration+Strata+Habitat,lambda="ML")

fPC1modsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  fPC1modsAIC[i,2]<-fPC1mods[[i]]$aicc
  fPC1modsAIC[i,3]<-summary(fPC1mods[[i]])$adj.r.squared
  fPC1modsAIC[i,4]<-pf(summary(fPC1mods[[i]])$fstatistic[1],summary(fPC1mods[[i]])$fstatistic[2],summary(fPC1mods[[i]])$fstatistic[3],lower.tail=F)
  fPC1modsAIC[i,5]<-fPC1mods[[i]]$param.CI$lambda$opt
  fPC1modsAIC[i,1]<-paste("fPC1_Model",i,sep="_")
}

fPC1modsAIC<-data.frame(fPC1modsAIC)
fPC1modsAIC[,2]<-as.numeric(as.character(fPC1modsAIC[,2]))
fPC1modsAIC[,3]<-as.numeric(as.character(fPC1modsAIC[,3]))
fPC1modsAIC[,5]<-as.numeric(as.character(fPC1modsAIC[,5]))

names(fPC1modsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda")

fPC1<-fPC1modsAIC[order(fPC1modsAIC$AICc),]
write.csv(fPC1,"Results/PGLS.nonphylo/fPC1_models.csv")

##############################
# Just categorical variables 
# Just categorical variables 
fPC2mods<-list()
fPC2mods[[1]]<-pgls(data=comp.data,fPC2~nonforest,lambda="ML")
fPC2mods[[2]]<-pgls(data=comp.data,fPC2~Strata,lambda="ML")
fPC2mods[[3]]<-pgls(data=comp.data,fPC2~Migration,lambda="ML")
fPC2mods[[4]]<-pgls(data=comp.data,fPC2~Forest.Dependency,lambda="ML")
fPC2mods[[5]]<-pgls(data=comp.data,fPC2~nonforest+Strata,lambda="ML")
fPC2mods[[6]]<-pgls(data=comp.data,fPC2~nonforest+Migration,lambda="ML")
fPC2mods[[7]]<-pgls(data=comp.data,fPC2~Forest.Dependency+Strata,lambda="ML")
fPC2mods[[8]]<-pgls(data=comp.data,fPC2~Forest.Dependency+Migration,lambda="ML")
fPC2mods[[9]]<-pgls(data=comp.data,fPC2~Migration+Strata,lambda="ML")
fPC2mods[[10]]<-pgls(data=comp.data,fPC2~Migration+Strata+nonforest,lambda="ML")
fPC2mods[[11]]<-pgls(data=comp.data,fPC2~Migration+Strata+Forest.Dependency,lambda="ML")
fPC2mods[[12]]<-pgls(data=comp.data,fPC2~Habitat,lambda="ML")
fPC2mods[[13]]<-pgls(data=comp.data,fPC2~Habitat+Strata,lambda="ML")
fPC2mods[[14]]<-pgls(data=comp.data,fPC2~Habitat+Migration,lambda="ML")
fPC2mods[[15]]<-pgls(data=comp.data,fPC2~Migration+Strata+Habitat,lambda="ML")

fPC2modsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  fPC2modsAIC[i,2]<-fPC2mods[[i]]$aicc
  fPC2modsAIC[i,3]<-summary(fPC2mods[[i]])$adj.r.squared
  fPC2modsAIC[i,4]<-pf(summary(fPC2mods[[i]])$fstatistic[1],summary(fPC2mods[[i]])$fstatistic[2],summary(fPC2mods[[i]])$fstatistic[3],lower.tail=F)
  fPC2modsAIC[i,5]<-fPC2mods[[i]]$param.CI$lambda$opt
  fPC2modsAIC[i,1]<-paste("fPC2_Model",i,sep="_")
}

fPC2modsAIC<-data.frame(fPC2modsAIC)
fPC2modsAIC[,2]<-as.numeric(as.character(fPC2modsAIC[,2]))
fPC2modsAIC[,3]<-as.numeric(as.character(fPC2modsAIC[,3]))
fPC2modsAIC[,5]<-as.numeric(as.character(fPC2modsAIC[,5]))

names(fPC2modsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda")

fPC2<-fPC2modsAIC[order(fPC2modsAIC$AICc),]
write.csv(fPC2,"Results/PGLS.nonphylo/fPC2_models.csv")

###########################################################################
## Brillance #########
# fAvgBrill
fAvgBrillmods<-list()
fAvgBrillmods[[1]]<-pgls(data=comp.data,fAvgBrill~nonforest,lambda="ML")
fAvgBrillmods[[2]]<-pgls(data=comp.data,fAvgBrill~Strata,lambda="ML")
fAvgBrillmods[[3]]<-pgls(data=comp.data,fAvgBrill~Migration,lambda="ML")
fAvgBrillmods[[4]]<-pgls(data=comp.data,fAvgBrill~Forest.Dependency,lambda="ML")
fAvgBrillmods[[5]]<-pgls(data=comp.data,fAvgBrill~nonforest+Strata,lambda="ML")
fAvgBrillmods[[6]]<-pgls(data=comp.data,fAvgBrill~nonforest+Migration,lambda="ML")
fAvgBrillmods[[7]]<-pgls(data=comp.data,fAvgBrill~Forest.Dependency+Strata,lambda="ML")
fAvgBrillmods[[8]]<-pgls(data=comp.data,fAvgBrill~Forest.Dependency+Migration,lambda="ML")
fAvgBrillmods[[9]]<-pgls(data=comp.data,fAvgBrill~Migration+Strata,lambda="ML")
fAvgBrillmods[[10]]<-pgls(data=comp.data,fAvgBrill~Migration+Strata+nonforest,lambda="ML")
fAvgBrillmods[[11]]<-pgls(data=comp.data,fAvgBrill~Migration+Strata+Forest.Dependency,lambda="ML")
fAvgBrillmods[[12]]<-pgls(data=comp.data,fAvgBrill~Habitat,lambda="ML")
fAvgBrillmods[[13]]<-pgls(data=comp.data,fAvgBrill~Habitat+Strata,lambda="ML")
fAvgBrillmods[[14]]<-pgls(data=comp.data,fAvgBrill~Habitat+Migration,lambda="ML")
fAvgBrillmods[[15]]<-pgls(data=comp.data,fAvgBrill~Migration+Strata+Habitat,lambda="ML")

fAvgBrillmodsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  fAvgBrillmodsAIC[i,2]<-fAvgBrillmods[[i]]$aicc
  fAvgBrillmodsAIC[i,3]<-summary(fAvgBrillmods[[i]])$adj.r.squared
  fAvgBrillmodsAIC[i,4]<-pf(summary(fAvgBrillmods[[i]])$fstatistic[1],summary(fAvgBrillmods[[i]])$fstatistic[2],summary(fAvgBrillmods[[i]])$fstatistic[3],lower.tail=F)
  fAvgBrillmodsAIC[i,5]<-fAvgBrillmods[[i]]$param.CI$lambda$opt
  #fAvgBrillmodsAIC[i,6]<-as.character(fAvgBrillmods[[i]]$varNames)
  fAvgBrillmodsAIC[i,1]<-paste("fAvgBrill_Model",i,sep="_")
}

fAvgBrillmodsAIC<-data.frame(fAvgBrillmodsAIC)
fAvgBrillmodsAIC[,2]<-as.numeric(as.character(fAvgBrillmodsAIC[,2]))
fAvgBrillmodsAIC[,3]<-as.numeric(as.character(fAvgBrillmodsAIC[,3]))
fAvgBrillmodsAIC[,5]<-as.numeric(as.character(fAvgBrillmodsAIC[,5]))
#fAvgBrillmodsAIC[,6]<-(as.character(fAvgBrillmodsAIC[,6]))
names(fAvgBrillmodsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda","traits")

fAvgBrill<-fAvgBrillmodsAIC[order(fAvgBrillmodsAIC$AICc),]
write.csv(fAvgBrill,"Results/PGLS.nonphylo/fAvgBrill_models.csv")
#############

# mAvgBrill
mAvgBrillmods<-list()
mAvgBrillmods[[1]]<-pgls(data=comp.data,mAvgBrill~nonforest,lambda="ML")
mAvgBrillmods[[2]]<-pgls(data=comp.data,mAvgBrill~Strata,lambda="ML")
mAvgBrillmods[[3]]<-pgls(data=comp.data,mAvgBrill~Migration,lambda="ML")
mAvgBrillmods[[4]]<-pgls(data=comp.data,mAvgBrill~Forest.Dependency,lambda="ML")
mAvgBrillmods[[5]]<-pgls(data=comp.data,mAvgBrill~nonforest+Strata,lambda="ML")
mAvgBrillmods[[6]]<-pgls(data=comp.data,mAvgBrill~nonforest+Migration,lambda="ML")
mAvgBrillmods[[7]]<-pgls(data=comp.data,mAvgBrill~Forest.Dependency+Strata,lambda="ML")
mAvgBrillmods[[8]]<-pgls(data=comp.data,mAvgBrill~Forest.Dependency+Migration,lambda="ML")
mAvgBrillmods[[9]]<-pgls(data=comp.data,mAvgBrill~Migration+Strata,lambda="ML")
mAvgBrillmods[[10]]<-pgls(data=comp.data,mAvgBrill~Migration+Strata+nonforest,lambda="ML")
mAvgBrillmods[[11]]<-pgls(data=comp.data,mAvgBrill~Migration+Strata+Forest.Dependency,lambda="ML")
mAvgBrillmods[[12]]<-pgls(data=comp.data,mAvgBrill~Habitat,lambda="ML")
mAvgBrillmods[[13]]<-pgls(data=comp.data,mAvgBrill~Habitat+Strata,lambda="ML")
mAvgBrillmods[[14]]<-pgls(data=comp.data,mAvgBrill~Habitat+Migration,lambda="ML")
mAvgBrillmods[[15]]<-pgls(data=comp.data,mAvgBrill~Migration+Strata+Habitat,lambda="ML")

mAvgBrillmodsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  mAvgBrillmodsAIC[i,2]<-mAvgBrillmods[[i]]$aicc
  mAvgBrillmodsAIC[i,3]<-summary(mAvgBrillmods[[i]])$adj.r.squared
  mAvgBrillmodsAIC[i,4]<-pf(summary(mAvgBrillmods[[i]])$fstatistic[1],summary(mAvgBrillmods[[i]])$fstatistic[2],summary(mAvgBrillmods[[i]])$fstatistic[3],lower.tail=F)
  mAvgBrillmodsAIC[i,5]<-mAvgBrillmods[[i]]$param.CI$lambda$opt
  mAvgBrillmodsAIC[i,1]<-paste("mAvgBrill_Model",i,sep="_")
}

mAvgBrillmodsAIC<-data.frame(mAvgBrillmodsAIC)
mAvgBrillmodsAIC[,2]<-as.numeric(as.character(mAvgBrillmodsAIC[,2]))
mAvgBrillmodsAIC[,3]<-as.numeric(as.character(mAvgBrillmodsAIC[,3]))
mAvgBrillmodsAIC[,5]<-as.numeric(as.character(mAvgBrillmodsAIC[,5]))

names(mAvgBrillmodsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda")

mAvgBrill<-mAvgBrillmodsAIC[order(mAvgBrillmodsAIC$AICc),]
write.csv(mAvgBrill,"Results/PGLS.nonphylo/mAvgBrill_models.csv")

##############################
######## Raw Traits ########
#############################
WPTCS <- read.csv("Colorspace/WPTCS.csv", header = TRUE)
female_WPTCS <- subset(WPTCS, WPTCS$Sex == "Female")
male_WPTCS <- subset(WPTCS, WPTCS$Sex == "Male")

rownames(male_WPTCS) <- male_WPTCS[,1,]
male_WPTCS<- male_WPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
male_WPTCS<-male_WPTCS[Card.tree$tip.label, ]  
male_WPTCS <- male_WPTCS %>% 
  rename(
    mAvgSpan = AvgSpan,
    mVarSpan = VarSpan,
    mMaxSpan = MaxSpan,
    mVolume = Volume,
    mAvgHueDisp = AvgHueDisp,
    mVarHueDisp = VarHueDisp,
    mMaxHueDisp = MaxHueDisp,
    mAvgBrill = AvgBrill,
    mAvgChroma = AvgChroma,
    mAvgAchChroma = AvgAchChroma)

rownames(female_WPTCS) <- female_WPTCS[,1,]
female_WPTCS<- female_WPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
female_WPTCS<-female_WPTCS[Card.tree$tip.label, ]

female_WPTCS <- female_WPTCS %>% 
  rename(
    fAvgSpan = AvgSpan,
    fVarSpan = VarSpan,
    fMaxSpan = MaxSpan,
    fVolume = Volume,
    fAvgHueDisp = AvgHueDisp,
    fVarHueDisp = VarHueDisp,
    fMaxHueDisp = MaxHueDisp,
    fAvgBrill = AvgBrill,
    fAvgChroma = AvgChroma,
    fAvgAchChroma = AvgAchChroma)


log.m_WPTCS <- log(male_WPTCS)
log.f_WPTCS <- log(female_WPTCS)



Color <- cbind(male_WPTCS,female_WPTCS)
Color$phylo <- row.names(Color)
comp.data<-comparative.data(Card.tree, Color, names.col = "phylo", vcv=TRUE, warn.dropped=TRUE)


FvMmods <-list()
FvMmods[[1]]<-pgls(data=comp.data,fAvgSpan~mAvgSpan,lambda="ML")
FvMmods[[2]]<-pgls(data=comp.data,fVarSpan~mVarSpan,lambda="ML")
FvMmods[[3]]<-pgls(data=comp.data,fMaxSpan~mMaxSpan,lambda="ML")
FvMmods[[4]]<-pgls(data=comp.data,fVolume~mVolume,lambda="ML")
FvMmods[[5]]<-pgls(data=comp.data,fAvgHueDisp~mAvgHueDisp,lambda="ML")
FvMmods[[6]]<-pgls(data=comp.data,fVarHueDisp~mVarHueDisp,lambda="ML")
FvMmods[[7]]<-pgls(data=comp.data,fMaxHueDisp~mMaxHueDisp,lambda="ML")
FvMmods[[8]]<-pgls(data=comp.data,fAvgBrill~mAvgBrill,lambda="ML")
FvMmods[[9]]<-pgls(data=comp.data,fAvgChroma~mAvgChroma,lambda="ML")
FvMmods[[10]]<-pgls(data=comp.data,fAvgAchChroma~mAvgAchChroma,lambda="ML")

FvMmodsAIC<-matrix(nrow=10,ncol=5)
for (i in 1:10){
  FvMmodsAIC[i,2]<-FvMmods[[i]]$aicc
  FvMmodsAIC[i,3]<-summary(FvMmods[[i]])$adj.r.squared
  FvMmodsAIC[i,4]<-pf(summary(FvMmods[[i]])$fstatistic[1],summary(FvMmods[[i]])$fstatistic[2],summary(FvMmods[[i]])$fstatistic[3],lower.tail=F)
  FvMmodsAIC[i,5]<-FvMmods[[i]]$param.CI$lambda$opt
  FvMmodsAIC[i,3]<-summary(FvMmods[[i]])$adj.r.squared
  FvMmodsAIC[i,1]<-paste("FvM_model",i,sep="_")
}
summary(FvMmods[[i]]$Intercept)

FvMmodsAIC<-data.frame(FvMmodsAIC)
FvMmodsAIC[,2]<-as.numeric(as.character(FvMmodsAIC[,2]))
FvMmodsAIC[,3]<-as.numeric(as.character(FvMmodsAIC[,3]))
FvMmodsAIC[,5]<-as.numeric(as.character(FvMmodsAIC[,5]))

names(FvMmodsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda")

FvM<-FvMmodsAIC[order(FvMmodsAIC$AICc),]
write.csv(FvM,"Results/PGLS.nonphylo/Chapter2/FvM_rawTraits_models.csv")


####################################################################################
####################################################################################
####################################################################################
####################################################################################

###################################################################
pdf("Results/PGLS.nonphylo/mPC2~VisRed.mean.pdf")


ggplot(comp.data$data,aes(x=mPC2,y=VisRed.mean))+geom_jitter(width=.1,height=.1,size=4,alpha=.5)+theme_bw()+geom_abline()

dev.off()


library(reshape)
d<-data.frame(comp.data$data$Percent.Tree.cover,comp.data$data$ClimPC1,comp.data$data$LandComp.1,comp.data$data$NDVI.mean)
names(d)<-c("Percent.Tree.cover","ClimPC1","LandComp.1","NDVI.mean")
td16<-melt(d,id="Percent.Tree.cover")

pdf("Results/PGLS.nonphylo/fPC1~Tree.Cover+Clim1+Land1.pdf")
ggplot(td16,aes(y=tdextent,x=value))+geom_jitter()+theme_bw()+geom_smooth(method="lm",alpha=.2,colour="gray")+facet_wrap(~variable,scales="free")
dev.off()
