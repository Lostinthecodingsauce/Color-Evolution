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

# Read in Tree and data 
Card.tree <- read.nexus("FullMCC.tree.nexus")
plot(Card.tree)

data <- read.csv("Chap2data.csv", header = TRUE) 

rownames(data) <- data[,1,]
Brill <- log(data[c(13,14)])
data <- data[c(-13,-14)]
df <- cbind(data,Brill)

model_names <- c("Nonforest","Strata","Migration","Forest Dependency","Nonforest+Strata","Nonforest+Migration",
                "Forest Dependency+Strata","Forest Dependency+Migration","Migration+Strata",
                 "Migration+Strata+Nonforest","Migration+Strata+Forest Dependency","Habitat",
                   "Habitat+Strata","Habitat+Migration","Migration+Strata+Habitat"
                   )

#Model_1	Nonforest
#Model_2	Strata
##Model_3	Migration
#Model_4	Forest Dependency
#Model_5	Nonforest+Strata
#Model_6	Nonforest+Migration
#Model_7	Forest Dependency+Strata
#Model_8	Forest Dependency+Migration
#Model_9	Migration+Strata
#Model_10	Migration+Strata+Nonforest
#Model_11	Migration+Strata+Forest Dependency
#Model_12	Habitat
#Model_13	Habitat+Strata
#Model_14	Habitat+Migration
#Model_15	Migration+Strata+Habitat


#Align species name to tree
df<-df[Card.tree$tip.label, ]
comp.data<-comparative.data(Card.tree, df, names.col = "phylo", vcv=TRUE, warn.dropped=TRUE)

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


###############   Female vs Male  ###############

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


###############   mPC1  ###############

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

mPC1modsAIC["Habitat_Trait"] <- c(model_names)
mPC1modsAIC["Plumage_Trait"] <- c("mPC1")

mPC1<-mPC1modsAIC[order(mPC1modsAIC$AICc),]


###############   mPC2  ###############

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
mPC2modsAIC["Habitat_Trait"] <- c(model_names)
mPC2modsAIC["Plumage_Trait"] <- c("mPC2")

mPC2<-mPC2modsAIC[order(mPC2modsAIC$AICc),]

###############   fPC1  ###############

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

fPC1modsAIC["Habitat_Trait"] <- c(model_names)
fPC1modsAIC["Plumage_Trait"] <- c("fPC1")

fPC1<-fPC1modsAIC[order(fPC1modsAIC$AICc),]

###############   fPC2  ############### 

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


fPC2modsAIC["Habitat_Trait"] <- c(model_names)
fPC2modsAIC["Plumage_Trait"] <- c("fPC2")

fPC2<-fPC2modsAIC[order(fPC2modsAIC$AICc),]

###############   F Brillance  ###############

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

fAvgBrillmodsAIC["Habitat_Trait"] <- c(model_names)
fAvgBrillmodsAIC["Plumage_Trait"] <- c("fBrill")

fAvgBrill<-fAvgBrillmodsAIC[order(fAvgBrillmodsAIC$AICc),]

###############   M Brillance  ###############

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

mAvgBrillmodsAIC["Habitat_Trait"] <- c(model_names)
mAvgBrillmodsAIC["Plumage_Trait"] <- c("mBrill")

mAvgBrill<-mAvgBrillmodsAIC[order(mAvgBrillmodsAIC$AICc),]

###################################
# Combine all results into single model

all_Merged <- do.call("rbind", list(mPC1, mPC2, fPC1, fPC2))#, mAvgBrill, fAvgBrill))

##################
