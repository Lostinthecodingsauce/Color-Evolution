##### Bayestraits vs ES-sims


library(ape)
library(mvtnorm)
library(geiger)
library(nlme)
library(phytools)
library(picante)
library(BSDA)
library(maps)
library(digest)
library(phylolm)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library (stats4)

setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis")
source("Code/RRPhylo.R")
source("Code/JetzDivRates.R")
source("Code/essim.R")
source("Code/matchplotdatasets.R") #This is a quick function I wrote to take two datasets and a phylogeny, and prune the datasets to include just the species they have in common, and plot them against each other 
source("Code/plotRates.R")


#Read in phenotypic data and tree------------
Cardinal.tree <- read.nexus("Cardinalidae.nexus")
FullCard.tree <- read.nexus("Cardinalidae.nexus")

All <- read.csv("All_reduced.csv")
m_traits <- read.csv("R_m_score.csv")
f_traits <- read.csv("R_f_score.csv")
mrates <- read.csv("MaleBayesTraits.csv")
frates <- read.csv("FemaleBayesTraits.csv")
rownames(mrates) <- mrates[,1,]
mrates<- mrates[,c(2,3)]

#tip <- c("Cyanocompsa_parellina")
#Cardinal.tree <- drop.tip(Cardinal.tree,tip)
#mrates<-mrates[Cardinal.tree$tip.label, ] 
#frates<-frates[Cardinal.tree$tip.label, ] 

###########################################################################################################################
# Male Traits
###########################################################################################################################
###################
#ES-sim on traits#
##################
traits_to_test <- c("PC1","PC2")
traits_m_res <- list()
traits_f_res <- list()
pdf("Color_Trait_Plots_M_PCA.pdf")

for (i in 1:length(traits_to_test)){
  #Male trait
  m_trait <- m_traits[,(traits_to_test[i])]
  names(m_trait) <- m_traits$Species 
  #mtrait <- log(m_trait)
  x<-treedata(Cardinal.tree, m_trait, warnings = FALSE)
  tree<-x$phy
  eqsp<-jetzDivRates(tree)
  eq.sp<-jetzDivRates(FullCard.tree)
  mtrait.p<-x$data
  traits_m_res[[i]] <- essim(trait = m_trait, phy = tree, nsim=1000, eq.sp)
  plot(mtrait.p, log(eqsp), main=paste0("male ",traits_to_test[[i]]))
}
dev.off()

##############
#RRPhyloRates#
##############
######### Testing whether tip-specific rates are associated with speciation rates ###########

#This code takes a trait, calculates tip-specific evolutionary rates with RRPhylo, and tests whether those rates are associated with speciation rates using ES-sim
#Because the rates produced by RRPhylo include both magnitude and direction (+/-), I used the absolute value of the rates
#A significant p-value means the absolute value of the rate (i.e. the magnitude) is significantly correlated with speciation rates

traits_to_test <- c("PC1","PC2")
traits_m_Bayes_res <- list() #Save the ES-sim results
rates_m_<-list() #Save the rates for tips
branchrates_m_rrphylo<-list() #save the rates for all edges in tree
pdf("BayesTraits_Tip_Rate_Plots_M_PCA.pdf") #Save plots of rates x equal splits 

for (i in 1:length(traits_to_test)){
  #Male trait
  m_trait <- m_traits[,(traits_to_test[i])]
  names(m_trait) <- m_traits$Species 
  #m_trait<-log(m_trait)
  m_trait<-m_trait[!is.na(m_trait)]
  #Match the tree and data
  x<-treedata(Cardinal.tree, m_trait, warnings=FALSE)
  tree<-x$phy
  data<-x$data
  #Calculate rates
  rates<-mrates
  #Extract all rates
  #allrates<-RRrates_trait$rates
 # branchrates_m_rrphylo[[i]]<-allrates
  #Get just the rates for the tips
  #y<-treedata(Cardinal.tree, allrates, warnings=FALSE)
  #tree2<-y$phy
  #rates<-y$data[,1]
  #rates_m_rrphylo[[i]]<-rates #Save rates
  eq.sp<-jetzDivRates(FullCard.tree) # calculates DR across whole tree - more accurate than just calculating for pruned tree
  traits_m_Bayes_res[[i]]<-essim(tree, log(abs(rates)), nsim = 1000, is = eq.sp)
  eqsp<-jetzDivRates(tree)
  plot<-plot(log(abs(rates)), log(eqsp), main=paste0("male ",traits_to_test[[i]]))
}

dev.off()

# Test for correlation
cor(eqsp, mrates$mPC1trate, method = "pearson")
