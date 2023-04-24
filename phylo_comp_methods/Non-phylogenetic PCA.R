# Non-phylogenetic PCA

library("phylopath")
library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("caper")
library(tidyverse)
library(ggplot2)
library(factoextra)
library(dplyr)
library(tidyverse)
library(glue)

setwd("C:/Research/Masters_thesis/MT_tree_analysis")


WPTCS <- read.csv("WPTCS-MT.csv", header = TRUE)
rownames(WPTCS) <- WPTCS[,1,]
WPTCS<- WPTCS[c(-1)]
log.WPTCS <- log(WPTCS[,1:9])
# Remove "Cyanoloxia_cyanoides")
#log.WPTCS.c <- log.WPTCS[-c(21), ] 
# Remove  "Piranga_hepatica"
#log.WPTCS.c <- log.WPTCS[-c(21,22,79,80,85,86), ] 

## Test for correlation amongst variables with chi-squared
chisq_results <- chisq.test() 

pca <- princomp(log.WPTCS, cor = TRUE, scores = TRUE)
pca.c <- princomp(log.WPTCS.c, cor = TRUE, scores = TRUE)
summary(pca)#Importance of Components
print(pca) #Standard Deviations

## Get scores
scores <- as.data.frame(pca$scores)
pca$loadings
pca.c$loadings



All.scores <- scores %>% 
  rename(
    AllPC.1 = Comp.1,
    AllPC.2 = Comp.2,
    AllPC.3 = Comp.3)


All.scores <- All.scores[c(1:3)]
All.scores$species <- row.names(All.scores)
sp<-ggplot(All.scores, aes(x=AllPC.1, y=AllPC.2, color=AllPC.3)) + geom_text(label=rownames(All.scores))
sp+scale_color_gradient(low="red", high="darkgreen")+ theme_minimal() 


#write.csv(All.scores, "raw_non_phylo.csv")
plotting <- read.csv("raw_non_phylo.csv")

plotting$species <- row.names(plotting)

     
ggplot(plotting) +
  aes(x = AllPC.1, y = AllPC.2, colour = sex) +
  geom_point(shape = "circle", size = 1.8) +
  scale_color_manual(values = c(F = "#6F98EE", M = "#B90F76")) +
  labs(
    x = "PC1 (56.1%)",
    y = "PC2 (28.5%)",
    title = "Nonphylogenetic PCA"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26L,
                              face = "bold",
                              hjust = 0.5),
    axis.title.y = element_text(size = 16L),
    axis.title.x = element_text(size = 16L)
  )



#write.csv(All.scores, file = "NonphyloScores.csv")
#########################################################
############### Plotting #########
fviz_eig(pca)
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )    # Avoid text overlapping
biplot.phyl.pca()
biplot(princomp(log.WPTCS, cor = TRUE, scores = TRUE))

var <- get_pca_var(pca)
var
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)
# Coordinates of variables
head(var$coord, 4)

fviz_pca_var(pca, col.var = "blue")

fviz_pca_var(pca, labelsize = 5, repel = TRUE) 

fviz_pca_var(pca, labelsize = 5, repel = TRUE)+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=96) +
  theme_minimal()



#######################################################################
####################    Begin PGLS     #################################
#################################################################
setwd("C:/Research/Masters_thesis/MT_tree_analysis")
#### PGLS for-loop #### 
library(ape)
library(phytools)
library(Rcpp)
library(geiger)
library(caper)
library(MCMCglmm)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(MuMIn)
  
# prepare phylogenetic datasets
#Card.tree <- read.nexus("Cardinalidae.nexus") #without Cyanoloxia_cyanoides
Card.tree <- read.nexus("FullMCC.tree.nexus")


df <- read.csv("Chap2data.nonphylo.csv", header = TRUE)
row.names(df) = df$Species

#get rid of node labels which prevent comparative dataframe object 
Card.tree$node.label<-NULL

comp.data<-comparative.data(Card.tree, df, names.col = "phylo", vcv=TRUE,vcv.dim = 3, warn.dropped=TRUE)
full <- comp.data$data

#full2 <- comp.data2$data
##############################################3
# sET UP LISTS FOR LATER
TRAITS <- c("Habitat", "Strata", "Migration","Forest Dependency","Openness+Strata","Habitat+Migration",
                                    "Forest Dependency+Strata","Forest Dependency+Migration","Migration+Strata",
                                      "Migration+Strata+Habitat","Migration+Strata+Forest Dependency","Openness",
                                        "Openness+Strata","Openness+Migration","Migration+Strata+Openness")
fpc1 < "fPC1"

###################
# Plots 
# Comparing male and female plumage 
df <- comp.data$data
a <- pgls(data=comp.data,mAvgBrill~Strata,lambda="ML")
summary(a)
a$aicc

with(df, plot(df, fPC1~mPC1, xlab = "fPC1", ylab = "mPC1",
                main = "fPC1 vs \ mPC1"))


sum <- summary(a)
sum$coefficients[2]
sum$coefficients[4]
#################
FvMmods<-list()
FvMmods[[1]]<-pgls(data=comp.data,fPC1~mPC1,lambda="ML")
FvMmods[[2]]<-pgls(data=comp.data,fPC2~mPC2,lambda="ML")
FvMmods[[3]]<-pgls(data=comp.data,fPC3~mPC3,lambda="ML")

FvMmodsAIC<-matrix(nrow=3,ncol=7)
for (i in 1:3){
  FvMmodsAIC[i,2]<-FvMmods[[i]]$aicc
  FvMmodsAIC[i,3]<-summary(FvMmods[[i]])$adj.r.squared
  FvMmodsAIC[i,4]<-pf(summary(FvMmods[[i]])$fstatistic[1],summary(FvMmods[[i]])$fstatistic[2],summary(FvMmods[[i]])$fstatistic[3],lower.tail=F)
  FvMmodsAIC[i,5]<-FvMmods[[i]]$param.CI$lambda$opt
  FvMmodsAIC[i,6]<-summary(FvMmods[[i]])$coefficients[2]
  FvMmodsAIC[i,7]<-summary(FvMmods[[i]])$coefficients[4]
  FvMmodsAIC[i,1]<-paste("FvM_Model",i,sep="_")
}

FvMmodsAIC<-data.frame(FvMmodsAIC)
FvMmodsAIC[,2]<-as.numeric(as.character(FvMmodsAIC[,2]))
FvMmodsAIC[,3]<-as.numeric(as.character(FvMmodsAIC[,3]))
FvMmodsAIC[,5]<-as.numeric(as.character(FvMmodsAIC[,5]))
FvMmodsAIC[,6]<-as.numeric(as.character(FvMmodsAIC[,6]))
FvMmodsAIC[,7]<-as.numeric(as.character(FvMmodsAIC[,7]))

names(FvMmodsAIC)<-c("model","AICc","adj_r_sq","p-value","lambda", "Beta", "Std_Error")

FvM<-FvMmodsAIC[order(FvMmodsAIC$AICc),]
write.csv(FvM,"C:/Research/Masters_thesis/MT_tree_analysis/Results/PGLS/Chapter2/nonphylo_pca/FvM_models.csv")


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

# Collect summary statistics
mPC1modsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  mPC1modsAIC[i,2]<-mPC1mods[[i]]$aicc
  mPC1modsAIC[i,3]<-summary(mPC1mods[[i]])$adj.r.squared
  mPC1modsAIC[i,4]<-pf(summary(mPC1mods[[i]])$fstatistic[1],summary(mPC1mods[[i]])$fstatistic[2],summary(mPC1mods[[i]])$fstatistic[3],lower.tail=F)
  mPC1modsAIC[i,5]<-mPC1mods[[i]]$param.CI$lambda$opt
  mPC1modsAIC[i,1]<-paste("Model",i,sep="_")
}

# transform to character
mPC1modsAIC<-data.frame(mPC1modsAIC)
mPC1modsAIC[,2]<-as.numeric(as.character(mPC1modsAIC[,2]))
mPC1modsAIC[,3]<-as.numeric(as.character(mPC1modsAIC[,3]))
mPC1modsAIC[,5]<-as.numeric(as.character(mPC1modsAIC[,5]))

# rename columns 
names(mPC1modsAIC)<-c("model","AICc","adj_r_sq","p_value","lambda")

# add plumage trait and habitat traits
mPC1modsAIC["Plumage_Trait"] <- c("mPC1")  ### change this out 
mPC1modsAIC["Habitat_Trait"] <- c(TRAITS)

mPC1modsAIC <- transform(mPC1modsAIC, `p_value` = as.numeric(`p_value`)) # convert p-value to numeric
mPC1modsAIC <- mPC1modsAIC %>% mutate(across(where(is.numeric), round, 3)) # round values to three decimal places

col_order <- c("Plumage_Trait", "Habitat_Trait", "model", "AICc", "adj_r_sq", "p_value", "lambda")
mPC1 <- mPC1modsAIC[, col_order]

mPC1<-mPC1[order(mPC1$AICc),]

write.csv(mPC1,"C:/Research/Masters_thesis/MT_tree_analysis/Results/PGLS/Chapter2/nonphylo_pca/mPC1_models.csv")

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

# Collect summary statistics
mPC2modsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  mPC2modsAIC[i,2]<-mPC2mods[[i]]$aicc
  mPC2modsAIC[i,3]<-summary(mPC2mods[[i]])$adj.r.squared
  mPC2modsAIC[i,4]<-pf(summary(mPC2mods[[i]])$fstatistic[1],summary(mPC2mods[[i]])$fstatistic[2],summary(mPC2mods[[i]])$fstatistic[3],lower.tail=F)
  mPC2modsAIC[i,5]<-mPC2mods[[i]]$param.CI$lambda$opt
  mPC2modsAIC[i,1]<-paste("Model",i,sep="_")
}

# transform to character
mPC2modsAIC<-data.frame(mPC2modsAIC)
mPC2modsAIC[,2]<-as.numeric(as.character(mPC2modsAIC[,2]))
mPC2modsAIC[,3]<-as.numeric(as.character(mPC2modsAIC[,3]))
mPC2modsAIC[,5]<-as.numeric(as.character(mPC2modsAIC[,5]))

# rename columns 
names(mPC2modsAIC)<-c("model","AICc","adj_r_sq","p_value","lambda")

# add plumage trait and habitat traits
mPC2modsAIC["Plumage_Trait"] <- c("mPC2")  ### change this out 
mPC2modsAIC["Habitat_Trait"] <- c(TRAITS)

mPC2modsAIC <- transform(mPC2modsAIC, `p_value` = as.numeric(`p_value`)) # convert p-value to numeric
mPC2modsAIC <- mPC2modsAIC %>% mutate(across(where(is.numeric), round, 3)) # round values to three decimal places

col_order <- c("Plumage_Trait", "Habitat_Trait", "model", "AICc", "adj_r_sq", "p_value", "lambda")
mPC2 <- mPC2modsAIC[, col_order]

mPC2<-mPC2[order(mPC2$AICc),]

write.csv(mPC2,"C:/Research/Masters_thesis/MT_tree_analysis/Results/PGLS/Chapter2/nonphylo_pca/mPC2_models.csv")

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

# Collect summary statistics
fPC1modsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  fPC1modsAIC[i,2]<-fPC1mods[[i]]$aicc
  fPC1modsAIC[i,3]<-summary(fPC1mods[[i]])$adj.r.squared
  fPC1modsAIC[i,4]<-pf(summary(fPC1mods[[i]])$fstatistic[1],summary(fPC1mods[[i]])$fstatistic[2],summary(fPC1mods[[i]])$fstatistic[3],lower.tail=F)
  fPC1modsAIC[i,5]<-fPC1mods[[i]]$param.CI$lambda$opt
  fPC1modsAIC[i,1]<-paste("Model",i,sep="_")
}

# transform to character
fPC1modsAIC<-data.frame(fPC1modsAIC)
fPC1modsAIC[,2]<-as.numeric(as.character(fPC1modsAIC[,2]))
fPC1modsAIC[,3]<-as.numeric(as.character(fPC1modsAIC[,3]))
fPC1modsAIC[,5]<-as.numeric(as.character(fPC1modsAIC[,5]))

# rename columns 
names(fPC1modsAIC)<-c("model","AICc","adj_r_sq","p_value","lambda")

# add plumage trait and habitat traits
fPC1modsAIC["Plumage_Trait"] <- c("fPC1")  ### change this out 
fPC1modsAIC["Habitat_Trait"] <- c(TRAITS)

fPC1modsAIC <- transform(fPC1modsAIC, `p_value` = as.numeric(`p_value`)) # convert p-value to numeric
fPC1modsAIC <- fPC1modsAIC %>% mutate(across(where(is.numeric), round, 3)) # round values to three decimal places

col_order <- c("Plumage_Trait", "Habitat_Trait", "model", "AICc", "adj_r_sq", "p_value", "lambda")
fPC1 <- fPC1modsAIC[, col_order]

fPC1<-fPC1[order(fPC1$AICc),]

write.csv(fPC1,"C:/Research/Masters_thesis/MT_tree_analysis/Results/PGLS/Chapter2/nonphylo_pca/fPC1_models.csv")

##############################
# Just categorical variables 
##############################
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

# Collect summary statistics
fPC2modsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  fPC2modsAIC[i,2]<-fPC2mods[[i]]$aicc
  fPC2modsAIC[i,3]<-summary(fPC2mods[[i]])$adj.r.squared
  fPC2modsAIC[i,4]<-pf(summary(fPC2mods[[i]])$fstatistic[1],summary(fPC2mods[[i]])$fstatistic[2],summary(fPC2mods[[i]])$fstatistic[3],lower.tail=F)
  fPC2modsAIC[i,5]<-fPC2mods[[i]]$param.CI$lambda$opt
  fPC2modsAIC[i,1]<-paste("Model",i,sep="_")
}

# transform to character
fPC2modsAIC<-data.frame(fPC2modsAIC)
fPC2modsAIC[,2]<-as.numeric(as.character(fPC2modsAIC[,2]))
fPC2modsAIC[,3]<-as.numeric(as.character(fPC2modsAIC[,3]))
fPC2modsAIC[,5]<-as.numeric(as.character(fPC2modsAIC[,5]))

# rename columns 
names(fPC2modsAIC)<-c("model","AICc","adj_r_sq","p_value","lambda")

# add plumage trait and habitat traits
fPC2modsAIC["Plumage_Trait"] <- c("fPC2")  ### change this out 
fPC2modsAIC["Habitat_Trait"] <- c(TRAITS)

fPC2modsAIC <- transform(fPC2modsAIC, `p_value` = as.numeric(`p_value`)) # convert p-value to numeric
fPC2modsAIC <- fPC2modsAIC %>% mutate(across(where(is.numeric), round, 3)) # round values to three decimal places

col_order <- c("Plumage_Trait", "Habitat_Trait", "model", "AICc", "adj_r_sq", "p_value", "lambda")
fPC2 <- fPC2modsAIC[, col_order]

fPC2<-fPC2[order(fPC2$AICc),]

write.csv(fPC2,"C:/Research/Masters_thesis/MT_tree_analysis/Results/PGLS/Chapter2/nonphylo_pca/fPC2_models.csv")

###########################################################################
## Brillance #########
# fAvgBrill
##############################
# Just categorical variables 
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

# Collect summary statistics
fAvgBrillmodsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  fAvgBrillmodsAIC[i,2]<-fAvgBrillmods[[i]]$aicc
  fAvgBrillmodsAIC[i,3]<-summary(fAvgBrillmods[[i]])$adj.r.squared
  fAvgBrillmodsAIC[i,4]<-pf(summary(fAvgBrillmods[[i]])$fstatistic[1],summary(fAvgBrillmods[[i]])$fstatistic[2],summary(fAvgBrillmods[[i]])$fstatistic[3],lower.tail=F)
  fAvgBrillmodsAIC[i,5]<-fAvgBrillmods[[i]]$param.CI$lambda$opt
  fAvgBrillmodsAIC[i,1]<-paste("Model",i,sep="_")
}

# transform to character
fAvgBrillmodsAIC<-data.frame(fAvgBrillmodsAIC)
fAvgBrillmodsAIC[,2]<-as.numeric(as.character(fAvgBrillmodsAIC[,2]))
fAvgBrillmodsAIC[,3]<-as.numeric(as.character(fAvgBrillmodsAIC[,3]))
fAvgBrillmodsAIC[,5]<-as.numeric(as.character(fAvgBrillmodsAIC[,5]))

# rename columns 
names(fAvgBrillmodsAIC)<-c("model","AICc","adj_r_sq","p_value","lambda")

# add plumage trait and habitat traits
fAvgBrillmodsAIC["Plumage_Trait"] <- c("fAvgBrill")  ### change this out 
fAvgBrillmodsAIC["Habitat_Trait"] <- c(TRAITS)

fAvgBrillmodsAIC <- transform(fAvgBrillmodsAIC, `p_value` = as.numeric(`p_value`)) # convert p-value to numeric
fAvgBrillmodsAIC <- fAvgBrillmodsAIC %>% mutate(across(where(is.numeric), round, 3)) # round values to three decimal places

col_order <- c("Plumage_Trait", "Habitat_Trait", "model", "AICc", "adj_r_sq", "p_value", "lambda")
fAvgBrill <- fAvgBrillmodsAIC[, col_order]

fAvgBrill<-fAvgBrill[order(fAvgBrill$AICc),]

write.csv(fAvgBrill,"C:/Research/Masters_thesis/MT_tree_analysis/Results/PGLS/Chapter2/nonphylo_pca/fAvgBrill_models.csv")
#############

# mAvgBrill
##############################
# Just categorical variables 
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

# Collect summary statistics
mAvgBrillmodsAIC<-matrix(nrow=15,ncol=5)
for (i in 1:15){
  mAvgBrillmodsAIC[i,2]<-mAvgBrillmods[[i]]$aicc
  mAvgBrillmodsAIC[i,3]<-summary(mAvgBrillmods[[i]])$adj.r.squared
  mAvgBrillmodsAIC[i,4]<-pf(summary(mAvgBrillmods[[i]])$fstatistic[1],summary(mAvgBrillmods[[i]])$fstatistic[2],summary(mAvgBrillmods[[i]])$fstatistic[3],lower.tail=F)
  mAvgBrillmodsAIC[i,5]<-mAvgBrillmods[[i]]$param.CI$lambda$opt
  mAvgBrillmodsAIC[i,1]<-paste("Model",i,sep="_")
}

# transform to character
mAvgBrillmodsAIC<-data.frame(mAvgBrillmodsAIC)
mAvgBrillmodsAIC[,2]<-as.numeric(as.character(mAvgBrillmodsAIC[,2]))
mAvgBrillmodsAIC[,3]<-as.numeric(as.character(mAvgBrillmodsAIC[,3]))
mAvgBrillmodsAIC[,5]<-as.numeric(as.character(mAvgBrillmodsAIC[,5]))

# rename columns 
names(mAvgBrillmodsAIC)<-c("model","AICc","adj_r_sq","p_value","lambda")

# add plumage trait and habitat traits
mAvgBrillmodsAIC["Plumage_Trait"] <- c("mAvgBrill")  ### change this out 
mAvgBrillmodsAIC["Habitat_Trait"] <- c(TRAITS)

mAvgBrillmodsAIC <- transform(mAvgBrillmodsAIC, `p_value` = as.numeric(`p_value`)) # convert p-value to numeric
mAvgBrillmodsAIC <- mAvgBrillmodsAIC %>% mutate(across(where(is.numeric), round, 3)) # round values to three decimal places

col_order <- c("Plumage_Trait", "Habitat_Trait", "model", "AICc", "adj_r_sq", "p_value", "lambda")
mAvgBrill <- mAvgBrillmodsAIC[, col_order]

mAvgBrill<-mAvgBrill[order(mAvgBrill$AICc),]

write.csv(mAvgBrill,"C:/Research/Masters_thesis/MT_tree_analysis/Results/PGLS/Chapter2/nonphylo_pca/mAvgBrill_models.csv")

############### Taking the best models ###################

df_list <- list(mPC1, mPC2, fPC1, fPC2, mAvgBrill, fAvgBrill)

best <- merge(df_list, aggregate(score ~ id, df_list=Plumage_traits, head, 1), by="id") 

# tidyverse/ FP solution
library(dplyr)
library(purrr)
library(readr)

map_df(.x = seq(1:3),
       .f = function(x) bind_rows(df_list[x]) %>%
         mutate(id = row_number(),
                df = x)) %>%
  arrange(id) %>%
  #select(-id, -df) %>% # uncomment if you want to lose row num and source
  write_csv(file = "C:/Research/Masters_thesis/MT_tree_analysis/Results/PGLS/Chapter2/nonphylo_pca/best_models.csv")
