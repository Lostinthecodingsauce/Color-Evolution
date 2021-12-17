#### Phylogenetic PCA on WPSTC ######
setwd("C:/Users/bfsco/Desktop/Masters Research/Chapter 2/analysis/WPTSC/MT_tree_analysis")

library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("RPANDA")
library(ggplot2)
library(phylosignal)
library(factoextra)
library(dplyr)

###### Organization #######
# 1) nonphylogenetic PCA both sexes
# 2) Male pPCA
# 3) female pPCA
# 4) Sexual Dichromatism attempt
# 7) Reduced plumage patches analysis, using Allisons plumage regions as well as belly and Medium Wing Coverts 

WPTCS <- read.csv("Colorspace/WPTCS.csv", header = TRUE)

female_WPTCS <- subset(WPTCS, WPTCS$Sex == "Female")
male_WPTCS <- subset(WPTCS, WPTCS$Sex == "Male")


Cardinal.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- read.nexus("FullMCC.tree.nexus")
FullCardinal.tree <- read.nexus("FullMCC.tree.nexus")




################################################################################
### 3) Phylogenetic PCA MALES
################################################################################
#male_WPTCS <- read.csv("Male.log.WPTCS.csv", header = TRUE)
rownames(male_WPTCS) <- male_WPTCS[,1,]
male_WPTCS<- male_WPTCS[,c(3,4,5,6,7,8,9,10,11,12)]
male_WPTCS<-male_WPTCS[Card.tree$tip.label, ]  

## Transform data 
log.m_WPTCS <- log(male_WPTCS[,1:10])


# log-transformed 
MPCA <- phyl.pca(Card.tree, log.m_WPTCS, method = "lambda", mode = "corr")

male_score <- as.data.frame(MPCA$S) ## Get scores
mload <- as.data.frame(MPCA$L)## Get loadings
summary(MPCA)#Importance of Components
malePRINT<- print(MPCA) #Standard Deviations

mPCA <-male_score[c(1:3)]

mPCA <- mPCA %>% 
  rename(
    mPC1 = PC1,
    mPC2 = PC2,
    mPC3 = PC3)

### Plot 
sp<-ggplot(mPCA, aes(x=mPC1, y=mPC2, color=mPC3)) + geom_text(label=rownames(mPCA))
sp+scale_color_gradient(low="red", high="darkgreen")

## Create martix for each PC score
malePC1<-as.matrix(male_score)[,1]
malePC2<-as.matrix(male_score)[,2]


## Map ancestral trait reconstruction 
anTPC1<-contMap(Cardinal.tree,malePC1, fsize=c(0.6,1),outline=FALSE)
anTPC2<-contMap(Cardinal.tree,malePC2, fsize=c(0.6,1),outline=FALSE)

#From the output, the p-value > 0.05 implying that the distribution of 
#the data are not significantly different from normal distribution.
#In other words, we can assume the normality.
shapiro.test(female_WPTCS$AvgSpan)
shapiro.test(female_WPTCS$VarSpan)
shapiro.test(female_WPTCS$MaxSpan) #no
shapiro.test(female_WPTCS$Volume) #no 
shapiro.test(female_WPTCS$AvgHueDisp) #no
shapiro.test(female_WPTCS$VarHueDisp) #no
shapiro.test(female_WPTCS$MaxHueDisp) #no
shapiro.test(female_WPTCS$AvgBrill) #no
shapiro.test(female_WPTCS$AvgChroma) #no
shapiro.test(female_WPTCS$AvgAchChroma) #no

hist(female_WPTCS$AvgSpan)
hist(female_WPTCS$VarSpan) # skewed right
hist(female_WPTCS$MaxSpan)
hist(female_WPTCS$Volume)
hist(female_WPTCS$AvgHueDisp)
hist(female_WPTCS$VarHueDisp)
hist(female_WPTCS$MaxHueDisp)
hist(female_WPTCS$AvgBrill)
hist(female_WPTCS$AvgChroma)
hist(female_WPTCS$AvgAchChroma)


################################################################################
### 3) Phylogenetic PCA FEMALES
################################################################################

rownames(female_WPTCS) <- female_WPTCS[,1,]
female_WPTCS<- female_WPTCS[,c(3,4,5,6,7,8,9,10,11,12)]

#Align species name to tree
female_WPTCS<-female_WPTCS[Card.tree$tip.label, ]

#log transform data
log.f_WPTCS <- log(female_WPTCS[,1:10])
 
FPCA <- phyl.pca(Card.tree, log.f_WPTCS, method = "lambda", mode = "corr")

summary(FPCA)#Importance of Components
print(FPCA) #Standard Deviations

## Get scores
female_score.log <- as.data.frame(FPCA$S)
female_load <- as.data.frame(FPCA$L)

fPCA <-female_score.log[c(1:3)]

fPCA <- fPCA %>% 
  rename(
    fPC1 = PC1,
    fPC2 = PC2,
    fPC3 = PC3)


PCscores <- cbind(mPCA,fPCA)
######################################################################

data <- cbind.data.frame(mPCA,fPCA)
data <- data*-1
data$phylo <- rownames(data)
write.csv(data,file = "Color.PC.scores.csv")

comp.data<-comparative.data(Card.tree, data, names.col = "phylo", vcv=TRUE, warn.dropped=TRUE)

a <- pgls(data=comp.data,fPC1~mPC1,lambda="ML")
summary(a)
with(data, plot(fPC1~mPC1, xlab = "mPC1", ylab = "fPC1",
                main = "Female PC1 vs \ Male PC1"))
abline(a)
# "Periporphyrus_erythromelas" tests as outlier

b <- pgls(data=comp.data,fPC2~mPC2,lambda="ML")
summary(b)
with(data, plot(fPC2~mPC2, xlab = "mPC2", ylab = "fPC2",
                main = "Female PC2 vs \ Male PC2"))
abline(b)

c <- pgls(data=comp.data,fPC3~mPC3,lambda="ML")
summary(c)
with(data, plot(fPC3~mPC3, xlab = "mPC3", ylab = "fPC3",
                main = "Female PC3 vs \ Male PC3"))
abline(c)

####################################################################
## Create martix for each PC score
femalePC1<-as.matrix(female_score)[,1]
femalePC2<-as.matrix(female_score)[,2]


## Map ancestral trait reconstruction 
anTPC1<-contMap(Cardinal.tree,femalePC1, fsize=c(0.6,1),outline=FALSE)
anTPC2<-contMap(Cardinal.tree,femalePC2, fsize=c(0.6,1),outline=FALSE)



## Map ancestral trait reconstruction 
anTPC1<-contMap(Card.tree,malePC1, fsize=c(0.6,1),outline=FALSE)
anTPC2<-contMap(Card.tree,malePC2, fsize=c(0.6,1),outline=FALSE)


#########################################

## Compare variance between raw traits of males and females using 2-sided T-test

AvgSpan <- var.test(male_WPTCS$AvgSpan, female_WPTCS$AvgSpan, alternative = "two.sided")
MaxSpan <- var.test(male_WPTCS$MaxSpan, female_WPTCS$MaxSpan, alternative = "two.sided")
VarSpan <- var.test(male_WPTCS$VarSpan, female_WPTCS$VarSpan, alternative = "two.sided")
Volume <- var.test(male_WPTCS$Volume, female_WPTCS$Volume, alternative = "two.sided")
AvgHue <- var.test(male_WPTCS$AvgHueDisp, female_WPTCS$AvgHueDisp, alternative = "two.sided")
VarHue<- var.test(male_WPTCS$VarHueDisp, female_WPTCS$VarHueDisp, alternative = "two.sided")
MaxHue<- var.test(male_WPTCS$MaxHueDisp, female_WPTCS$MaxHueDisp, alternative = "two.sided")
AvgBrill <- var.test(male_WPTCS$AvgBrill, female_WPTCS$AvgBrill, alternative = "two.sided")
AvgChroma<- var.test(male_WPTCS$AvgChroma, female_WPTCS$AvgChroma, alternative = "two.sided")
AvgAchChroma<- var.test(male_WPTCS$AvgAchChroma, female_WPTCS$AvgAchChroma, alternative = "two.sided")

t_test_WPTCS <- cbind(AvgAchChroma,AvgChroma,AvgBrill,MaxHue,VarHue,AvgHue,Volume,VarSpan,MaxSpan,AvgSpan)
t_test_WPTCS <-data.frame(t_test_WPTCS)

VarWPTCS <- var.test(male_WPTCS, female_WPTCS, alternative = "two.sided")
