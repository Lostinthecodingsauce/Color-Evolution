#Card ES-Sims_RRphylo.R

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
#source("Code/matchplotdatasets.R") #This is a quick function I wrote to take two datasets and a phylogeny, and prune the datasets to include just the species they have in common, and plot them against each other 
source("Code/plotRates.R")


#Read in phenotypic data and tree------------
tip <- c("Cyanocompsa_cyanoides")
#FullCard.tree <- read.nexus("MCC.Cardinalidae.nexus")
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Cardinal.tree <- drop.tip(FullCard.tree,tip)

all_color_data <- read.csv("WPTCS.csv", header = TRUE, na = "0") 
color_data_m <- read.csv("MT_Tree_Male_WPTCS.csv", header = TRUE)
color_data_f <- read.csv("MT_Tree_Female_WPTCS.csv", header = TRUE)
All <- read.csv("New.csv")
fuck <- read.csv("fuck.csv")

############################################# PCAs ############################################################3
#tip <- c("Cyanocompsa_cyanoides")
#Card.tree <- drop.tip(Cardinal.tree,tip)


###################
#ES-sim on traits#
##################
traits_to_test <- c("mPC1","mPC2")
traits_m_res <- list()
traits_f_res <- list()
pdf("mPCA_RRphylo.pdf")

for (i in 1:length(traits_to_test)){
  #Male trait
  m_trait <- All[,(traits_to_test[i])]
  names(m_trait) <- All$Species 
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

tip <- c("Cyanocompsa_cyanoides","Amaurospiza_concolor")
#FullCard.tree <- read.nexus("MCC.Cardinalidae.nexus")
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Cardinal.tree <- drop.tip(FullCard.tree,tip)

traits_to_test <- c("mPC1","mPC2")
traits_m_rrphylorates_res <- list() #Save the ES-sim results
rates_m_rrphylo<-list() #Save the rates for tips
branchrates_m_rrphylo<-list() #save the rates for all edges in tree
pdf("RRPhylo_Tip_Rate_Plots_M_PCA.pdf") #Save plots of rates x equal splits 

for (i in 1:length(traits_to_test)){
  #Male trait
  m_trait <- fuck[,(traits_to_test[i])]
  names(m_trait) <- fuck$ï..Species 
  #m_trait<-log(m_trait)
  m_trait<-m_trait[!is.na(m_trait)]
  #Match the tree and data
  x<-treedata(Cardinal.tree, m_trait, warnings=FALSE)
  tree<-x$phy
  data<-x$data
  #Calculate rates
  RRrates_trait<-RRphylo(tree, data)
  #Extract all rates
  allrates<-RRrates_trait$rates
  branchrates_m_rrphylo[[i]]<-allrates
  #Get just the rates for the tips
  y<-treedata(Cardinal.tree, allrates, warnings=FALSE)
  tree2<-y$phy
  rates<-y$data[,1]
  rates_m_rrphylo[[i]]<-rates #Save rates
  eq.sp<-jetzDivRates(FullCard.tree) # calculates DR across whole tree - more accurate than just calculating for pruned tree
  traits_m_rrphylorates_res[[i]]<-essim(tree, log(abs(rates)), nsim = 1000, is = eq.sp)
  eqsp<-jetzDivRates(tree)
  plot<-plot(log(abs(rates)), log(eqsp), main=paste0("male ",traits_to_test[[i]]))
}

dev.off()
##### Plot rates
a <- contMap(Cardinal.tree, rates_m_rrphylo)

### Female PCA#############
traits_to_test <- c("fPC1","fPC2")
traits_f_rrphylorates_res <- list() #Save the ES-sim results
rates_f_rrphylo<-list() #Save the rates for tips
branchrates_f_rrphylo<-list() #save the rates for all edges in tree
pdf("RRPhylo_Tip_Rate_Plots_F_PCA.pdf") #Save plots of rates x equal splits 

for (i in 1:length(traits_to_test)){
  #Female trait
  f_trait <- fuck[,(traits_to_test[i])]
  names(f_trait) <- fuck$ï..Species 
  #m_trait<-log(m_trait)
  f_trait<-f_trait[!is.na(f_trait)]
  #Match the tree and data
  x<-treedata(Cardinal.tree, f_trait, warnings=FALSE)
  tree<-x$phy
  data<-x$data
  #Calculate rates
  RRrates_trait<-RRphylo(tree, data)
  #Extract all rates
  allrates<-RRrates_trait$rates
  branchrates_f_rrphylo[[i]]<-allrates
  #Get just the rates for the tips
  y<-treedata(Cardinal.tree, allrates, warnings=FALSE)
  tree2<-y$phy
  rates<-y$data[,1]
  rates_f_rrphylo[[i]]<-rates #Save rates
  eq.sp<-jetzDivRates(FullCard.tree) # calculates DR across whole tree - more accurate than just calculating for pruned tree
  traits_f_rrphylorates_res[[i]]<-essim(tree, log(abs(rates)), nsim = 1000, is = eq.sp)
  eqsp<-jetzDivRates(tree)
  plot<-plot(log(abs(rates)), log(eqsp), main=paste0("female ",traits_to_test[[i]]))
}

dev.off()

########
woof <- read.csv("AAAAAAAA.csv", header = TRUE)
a <- contMap(Cardinal.tree, woof$mPC1)
