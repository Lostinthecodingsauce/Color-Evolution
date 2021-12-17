### Path analysuis 
library("phylopath")
library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("caper")
library(ggplot2)
# Examples https://cran.r-project.org/web/packages/phylopath/vignettes/intro_to_phylopath.html
# Examples https://peerj.com/articles/4718/

# Introduction to PPA
# http://www.mpcm-evolution.com/practice/online-practical-material-chapter-8/chapter-8-2-step-step-guide-phylogenetic-path-analysis-using-d-sep-method-rhinograds-example


setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis")
#tip <- c("Cyanocompsa_cyanoides","Caryothraustes_canadensis", "Amaurospiza_concolor")
#tip <- c("Cyanocompsa_cyanoides","Caryothraustes_canadensis")
tip <- c("Cyanocompsa_cyanoides")
#FullCard.tree <- read.nexus("MCC.Cardinalidae.nexus")
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- drop.tip(FullCard.tree,tip)

All <- read.csv("New.csv", header = TRUE) # without Cyanocompsa cyanoides
rownames(All) <- This[,1,]
All<- All[,c(2,3,4,6,7,8,9,10,11,12,13,14,15)]
All<-All[Card.tree$tip.label, ]  
#### Extract columns,
f_PC1 <- All[, "fPC1A"]
f_PC2 <- All[, "fPC2A"]
m_PC1 <- All[, "mPC1"]
m_PC2 <- All[, "mPC2"]
Habitat <- All[, "Habitat"]
Strata <- All[, "StrataA"]
Mig <- All[,"Migration"]
Fdep <- All[,"Forest.Dependency"]
Forest <- All[, "FUCK"]

#give them names
names(Forest) <- names(Mig) <- names(Strata) <- names(Habitat) <- names(m_PC2) <- names(m_PC1) <- rownames(All)
names(f_PC2) <- names(f_PC1)<- names(Fdep)  <- rownames(All)

models <- define_model_set(
  a = c(mPC1 ~ fPC1),
  b = c(mPC1 ~ fPC1+StrataA),
  c = c(mPC1 ~ fPC1+Habitat),
  d = c(mPC1 ~ fPC1+Mig),
  e = c(mPC1 ~ fPC1+Fdep),
  f = c(mPC1 ~ fPC1+Forest),
  g = c(mPC1 ~ fPC1 + StrataA+Habitat),
  h = c(mPC1 ~ fPC1 + StrataA+Forest),
  i = c(mPC1 ~ fPC1 + StrataA+Mig),
  j = c(mPC1 ~ fPC1 + StrataA+Fdeb),
  k = c(mPC1 ~ fPC1 + StrataA+Fdeb+Mig),
  l = c(mPC1 ~ fPC1 + StrataA+Fdeb+Habitat),
  m = c(mPC1 ~ fPC1 + StrataA+Fdeb+Forest),

  
  six   = c(NL ~ RS, RS ~ BM),
  seven = c(NL ~ RS, RS ~ LS + BM),
  eight = c(NL ~ RS),
  nine  = c(NL ~ RS, RS ~ LS),
  .common = c(LS ~ BM, NL ~ BM, DD ~ NL)
)



m <- define_model_set(
  null = c(),
  directMale = c(mPC1~fPC1),
  directFemale = c(fPC1~mPC1),
  indirectMale = c(Habitat~mPC1, StrataA~mPC1, Forest~mPC1, Migration~mPC1),
  indirectFemale = c(Habitat~fPC1, StrataA~fPC1, Forest~fPC1, Migration~fPC1),
  bothMale = c(mPC1~fPC1,Habitat~mPC1, StrataA~mPC1, Forest~mPC1, Migration~mPC1),
  bothFemale = c(fPC1~mPC1,Habitat~mPC1, StrataA~mPC1, Forest~mPC1, Migration~mPC1))
  #common = c(Br~B, P~B, L~B+G, W~G, Status~P+L+G+W+B))

plot_model_set(m)




m <- define_model_set(
  null = c(),
  direct = c(Status~Br),
  indirect = c(L~Br, G~Br, W~Br),
  both = c(Status~Br, L~Br, G~Br, W~Br),
  common = c(Br~B, P~B, L~B+G, W~G, Status~P+L+G+W+B))

plot_model_set(m)



### example 
models <- define_model_set(
  one   = c(RS ~ DD),
  two   = c(DD ~ NL, RS ~ LS + DD),
  three = c(RS ~ NL),
  four  = c(RS ~ BM + NL),
  five  = c(RS ~ BM + NL + DD),
  six   = c(NL ~ RS, RS ~ BM),
  seven = c(NL ~ RS, RS ~ LS + BM),
  eight = c(NL ~ RS),
  nine  = c(NL ~ RS, RS ~ LS),
  .common = c(LS ~ BM, NL ~ BM, DD ~ NL)
)
plot_model_set(models)
