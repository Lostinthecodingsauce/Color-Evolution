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

setwd("/Users/rosalynprice-waldman/Documents/Tanagers/PlumageSpeciationClean")
source("RRPhylo.R")
source("JetzDivRates.R")
source("essim.R")
source("matchplotdatasets.R") #This is a quick function I wrote to take two datasets and a phylogeny, and prune the datasets to include just the species they have in common, and plot them against each other 
source("plotRates.R")


#Read in phenotypic data and tree------------
tanagerTree <- read.tree("TanagerMCCclocked_edited.tre") #This tips on this tree should match the dataset

color_data <- read.csv("WPTCSmeasures.csv", header = TRUE, na = "0") 


###################
#ES-sim on traits#
##################
traits_to_test <- c("AvgSpan","Volume","AvgHueDisp","AvgBrill","AvgChroma")
traits_m_res <- list()
traits_f_res <- list()
pdf("Color_Trait_Plots.pdf")

for (i in 1:length(traits_to_test)){
  #Male trait
  m_trait <- color_data_m[,(traits_to_test[i])]
  names(m_trait) <- color_data_m$species
  mtrait <- log(m_trait)
  x<-treedata(tanagerTree, mtrait, warnings = FALSE)
  tree<-x$phy
  eqsp<-jetzDivRates(tree)
  eq.sp<-jetzDivRates(tanagerTree)
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

traits_to_test <- c("AvgSpan","Volume","AvgHueDisp","AvgBrill","AvgChroma")
traits_m_rrphylorates_res <- list() #Save the ES-sim results
rates_m_rrphylo<-list() #Save the rates for tips
branchrates_m_rrphylo<-list() #save the rates for all edges in tree
pdf("RRPhylo_Tip_Rate_Plots.pdf") #Save plots of rates x equal splits 

for (i in 1:length(traits_to_test)){
  #Male trait
  m_trait <- color_data_m[,(traits_to_test[i])]
  names(m_trait) <- color_data_m$species
  m_trait<-log(m_trait)
  m_trait<-m_trait[!is.na(m_trait)]
  #Match the tree and data
  x<-treedata(tanagerTree, m_trait, warnings=FALSE)
  tree<-x$phy
  data<-x$data
  #Calculate rates
  RRrates_trait<-RRphylo(tree, data)
  #Extract all rates
  allrates<-RRrates_trait$rates
  branchrates_m_rrphylo[[i]]<-allrates
  #Get just the rates for the tips
  y<-treedata(tanagerTree, allrates, warnings=FALSE)
  tree2<-y$phy
  rates<-y$data[,1]
  rates_m_rrphylo[[i]]<-rates #Save rates
  eq.sp<-jetzDivRates(tanagerTree) # calculates DR across whole tree - more accurate than just calculating for pruned tree
  traits_m_rrphylorates_res[[i]]<-essim(tree, log(abs(rates)), nsim = 1000, is = eq.sp)
  eqsp<-jetzDivRates(tree)
  plot(log(abs(rates)), log(eqsp), main=paste0("male ",traits_to_test[[i]]))
}

dev.off()