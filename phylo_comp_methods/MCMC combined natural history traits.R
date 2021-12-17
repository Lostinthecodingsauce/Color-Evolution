#MCMCglmm

install.packages("MCMCglmm")

library("brms")
library("MCMCglmm")


setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/analysis/WPTSC/MT_tree_analysis")

#FullCard.tree <- read.nexus("MCC.Cardinalidae.nexus")
Card.tree <- read.nexus("Cardinalidae.nexus")
A <- ape::vcv.phylo(Card.tree)
This <- read.csv("New.csv", header = TRUE) # without Cyanocompsa cyanoides


model_simple <- brm(
  mPC1 ~ fPC1A + (1|Card.tree), data = This, 
  family = gaussian(), cov_ranef = list(phylo = A),
  prior = c(
    prior(normal(0, 10), "b"),
    prior(normal(0, 50), "Intercept"),
    prior(student_t(3, 0, 20), "sd"),
    prior(student_t(3, 0, 20), "sigma")
  )
)


phylo <- ape::read.nexus("Cardinalidae.nexus")
data_simple <- read.table("New.csv", header = TRUE)
head(data_simple)