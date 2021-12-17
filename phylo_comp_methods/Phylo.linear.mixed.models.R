# Phylogenetic Linear Mixed Models 
## https://github.com/paul-buerkner/brms
## https://cran.r-project.org/web/packages/MCMCglmm/vignettes/CourseNotes.pdf
## https://github.com/tmalsburg/MCMCglmm-intro

#MCMCglmm walkthrough
# http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

library(ape)
library(phytools)
library(Rcpp)
library(geiger)
library(RPANDA)
library(caper)
library(MCMCglmm)
library(ggplot2)

########################################################################

setwd("C:/Users/bfsco/Desktop/Masters Research/Plumage/analysis/WPTSC/MT_tree_analysis")
#tip <- c("Cyanocompsa_cyanoides") 
FullCard.tree <- read.nexus("FullMCC.tree.nexus")#, force.multi = TRUE)
Card.tree <- read.nexus("FullMCC.tree.nexus")
#Card.tree <- drop.tip(FullCard.tree,tip)

# Multiple different data sets. This one has PCAs that were not log-transformed, all non-scaled spatial data, 
#data <- read.csv("MCMC_NDVI.new.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') # PCA not log-transformed 
newer <- read.csv("IUCN/new.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') # PCA not log-transformed 

# Transform Spatial data to have correct scale based on NDVI values 
df.MODIS <- newer[c(16:21)]*0.0001 # df.mean <- newer[c(16,21,26,31)]*0.0001
df.GPP <- newer[c(22:23)]*0.01 # df.mean <- newer[c(16,21,26,31)]*0.0001


temp <- subset(newer[c(1:15,24:26,31:32)])# subset(All[c(1:15,17:20,22:25,27:30,32:36)])
temp2 <- cbind.data.frame(temp,df.MODIS)
data <- cbind.data.frame(temp2,df.GPP)


#comp.data<-comparative.data(Card.tree, data, names.col = "phylo", vcv=TRUE, warn.dropped=TRUE)
comp.data<-comparative.data(Card.tree, data, names.col = "X", vcv=TRUE, warn.dropped=TRUE)


full <- data
rownames(full) <- full[,1,]
full<- full[-c(1)]
########################


################################################################################
################### Bayesian Analysis #################################

inv.phylo<-inverseA(Card.tree,nodes="TIPS",scale=TRUE)

# Define priors. Use default prior, test other priors later 
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

model_simple<-MCMCglmm(fPC1~mPC1,random=~phylo,
                       family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                       prior=prior,data=All,nitt=500000,burnin=1000,thin=500)

#Summary
summary(model_simple)
plot(model_simple)
lambda <- model_simple$VCV[,'phylo']/
  (model_simple$VCV[,'phylo']+model_simple$VCV[,'units'])
mean(lambda)
## [1] 0.6961
posterior.mode(lambda)
##   var1 
## 0.7442
HPDinterval(lambda)
############################### Works !!!!!!!!!!!!!!!!!!
# Testing using multiple traits 
#NDVI.mean	NDVI.stdDev	NDVI.min	NDVI.max	NDVI.skew	
#EVI.mean	EVI.stdDev	EVI.min	EVI.max	EVI.skew



model_m1_NDVImean<-MCMCglmm(mPC1~NDVI.mean,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_m1_NDVISD<-MCMCglmm(mPC1~NDVI.stdDev,random=~phylo,
                          family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                          prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)


model_f1_NDVImean<-MCMCglmm(fPC1~NDVI.mean,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_f1_NDVISD<-MCMCglmm(fPC1~NDVI.stdDev,random=~phylo,
                          family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                          prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)

####################################
model_m2_NDVImean<-MCMCglmm(mPC2~NDVI.mean,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_m2_NDVISD<-MCMCglmm(mPC2~NDVI.stdDev,random=~phylo,
                          family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                          prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)

#Female PC2
model_f2_NDVImean<-MCMCglmm(fPC2~NDVI.mean,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_f2_NDVISD<-MCMCglmm(fPC2~NDVI.stdDev,random=~phylo,
                          family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                          prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)


### Skew

model_m1_NDVISkew<-MCMCglmm(mPC1~NDVI.skew,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_m2_NDVISkew<-MCMCglmm(mPC2~NDVI.skew,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_f1_NDVISkew<-MCMCglmm(fPC1~NDVI.skew,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_f2_NDVISkew<-MCMCglmm(fPC2~NDVI.skew,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
#######################################################################

## Females
#NDVI
model_fB_NDVImean<-MCMCglmm(FemaleBrill~NDVI.mean,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_fB_NDVISD<-MCMCglmm(FemaleBrill~NDVI.stdDev,random=~phylo,
                          family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                          prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
#EVI
model_fB_EVImean<-MCMCglmm(FemaleBrill~EVI.mean,random=~phylo,
                           family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                           prior=prior,data=All,nitt=50000,burnin=10000,thin=1000)
model_fB_EVISD<-MCMCglmm(FemaleBrill~EVI.stdDev,random=~phylo,
                         family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                         prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
#NIR
model_fB_NIRmean<-MCMCglmm(FemaleBrill~NIR.mean,random=~phylo,
                           family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                           prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_fB_NIRSD<-MCMCglmm(FemaleBrill~NIR.stdDev,random=~phylo,
                         family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                         prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
#VisRed
model_fB_VisRedmean<-MCMCglmm(FemaleBrill~VisRed.mean,random=~phylo,
                              family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                              prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)
model_fB_VisRedSD<-MCMCglmm(FemaleBrill~VisRed.stdDev,random=~phylo,
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                            prior=prior,data=All,nitt=500000,burnin=10000,thin=1000)

###############
# Resuklts ##
###############


#Summary
summary(model_m1_NDVImean)
summary(model_f1_NDVImean)
summary(model_m2_NDVImean)
summary(model_f2_NDVImean)
# NDVI SD
summary(model_m1_NDVISD)
summary(model_f1_NDVISD)
summary(model_m2_NDVISD)
summary(model_f2_NDVISD)
# NDVI Skew
summary(model_m1_NDVISkew)
summary(model_f1_NDVISkew)
summary(model_m2_NDVISkew)
summary(model_f2_NDVISkew)
# EVI Mean f2_EVImean
summary(model_m1_EVImean)
summary(model_f1_EVImean)
summary(model_m2_EVImean)
summary(model_f2_EVImean)
# EVI SD 
summary(model_m1_EVISD)
summary(model_m2_EVISD)
summary(model_f1_EVISD)
summary(model_f2_EVISD)
# NIR
summary(model_m1_NIRmean)
summary(model_f1_NIRmean)
summary(model_m1_NIRSD)
summary(model_f1_NIRSD)
summary(model_m2_NIRmean)
summary(model_f2_NIRmean)
summary(model_m2_NIRSD)
summary(model_f2_NIRSD)


# Visable Red
summary(model_m1_VisRedmean)
summary(model_f1_VisRedmean)
summary(model_m1_VisRedSD)
summary(model_f1_VisRedSD)
summary(model_m2_VisRedmean)
summary(model_f2_VisRedmean)
summary(model_m2_VisRedSD)
summary(model_f2_VisRedSD)
##################################################
###### Brilliance 
summary(model_mB_NDVImean)
summary(model_fB_NDVImean)
summary(model_mB_NDVISD)
summary(model_fB_NDVISD)
summary(model_mB_EVImean)
summary(model_fB_EVImean) # signicant , but needs to be ran much longer 
summary(model_mB_EVISD)
summary(model_fB_EVISD)
summary(model_mB_NIRmean) # Significant 
summary(model_fB_NIRmean) # Significant 
summary(model_mB_NIRSD)
summary(model_fB_NIRSD)
summary(model_mB_VisRedmean) # Significant 
summary(model_fB_VisRedmean)# Significant 
summary(model_mB_VisRedSD) # Significant 
summary(model_fB_VisRedSD) # Significant 
################################################

####
plot(model_m1NDVImean)
lambda <- model_m1NDVImean$VCV[,'phylo']/
  (model_m1NDVImean$VCV[,'phylo']+model_m1NDVImean$VCV[,'units'])
mean(lambda)
## [1] 0.6961
posterior.mode(lambda)
##   var1 
## 0.7442
HPDinterval(lambda)


##############################################################################
###### Analysis of life history traits impact on plumage evolution


setwd("C:/Users/bfsco/Desktop/Masters Research/Plumage/analysis/WPTSC/MT_tree_analysis")
tip <- c("Cyanocompsa_cyanoides") #"Cyanocompsa_cyanoides"
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- drop.tip(FullCard.tree,tip)

All <- read.csv("Alldata.csv")
# Test 
# Construct covariance matrix of species
A <- ape::vcv.phylo(Card.tree)
# Built first model
library(lme4)
m1 <- glmer(pronoun ~  (a + b + c)^3            +
              ((a + b + c)^3 | subject) +
              ((a + b    )^2 | item),
            data=d, family="binomial")




model_1 <- brm(
  fPC1A ~ mPC1A + (1|gr(ï..species, cov = A)), 
  data = All, 
  family = gaussian(), 
  warmup = 2000, 
  iter   = 20000, 
  chains = 4, 
  data2 = list(A = A))

# get results
summary(model_1)

# PLot
plot(model_1, N = 2, ask = FALSE)
plot(conditional_effects(model_1), points = TRUE) 

### Check to see if it converged using catepillar plots
model1tranformed <- ggs(model1) # the ggs function transforms the brms output into a longformat tibble, that we can use to make different types of plots.
## Warning in custom.sort(D$Parameter): NAs introduced by coercion
ggplot(filter(model1tranformed, Parameter %in% c("b_Intercept", "b_extrav", "b_sex")),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = 1000)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Caterpillar Plots", 
       col   = "Chains")



##############################################################################
# Second model
model_2 <- brm(
  mPC1A ~ NDVI.stdDev + BioPC1 + (1|gr(ï..species, cov = A)), 
  data = This, 
  family = gaussian(), 
  #family = zero_inflated_poisson("log"),
  #family = lognormal(),
  data2 = list(A = A))

summary(model_2)
# PLot
plot(model_2, N = 2, ask = FALSE)
plot(conditional_effects(model_2), points = TRUE) 
marginal_effects(model_2)

# plotting?
library(ggstatsplot)
ggstatsplot::ggcoefstats(
  x = model_mixed,
  title = "MalePC1: Forest & NDVI.mean")


# PLot
plot(model_mixed, N = 2, ask = FALSE)
plot(conditional_effects(model_mixed), points = TRUE) 
marginal_effects(model_mixed)

# plotting?
library(ggstatsplot)
ggstatsplot::ggcoefstats(
  x = model_mixed,
  title = "MalePC2: Life History")
##############################################
# Plotting MCMC data
# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html#point-summaries-and-intervals
# https://www.rensvandeschoot.com/tutorials/brms-started/
# https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcoefstats.html

# https://rdrr.io/cran/phyr/man/pglmm.html
library(phyr)

#z <- pglmm(Y ~ X + (1|sp__), data = data, family = "gaussian", cov_ranef = list(sp = phy))
z <- pglmm(mPC1 ~ fPC1 + (1|sp__), data = All, family = "gaussian", 
           cov_ranef = list(sp = Card.tree))

modelflag <- 1
sim.dat <- PGLMM.sim(stree(16, "balanced"), nsites = 30,
                     modelflag = modelflag, + second.env = TRUE, compscale = 1)

#######################################################################
library(mvMORPH) #https://www.rdocumentation.org/packages/mvMORPH/versions/1.1.3
# Threshold model 
# Multivariate Analysis of Variance (manova.gls)

######################################################################
# Phylogenetic linear mixed models 

#https://cran.r-project.org/web/packages/pez/vignettes/pez-pglmm-overview.pdf

nspp <- 15
nsite <- 10
env <- 1:nsite
env <- as.numeric(scale(env))

#
require(pez)
require(ape)
phy <- rcoal(n=7)
Vphy <- vcv(phy)
Vphy <- Vphy/(det(Vphy)^(1/7))

iD <- t(chol(Vphy))
intercept <- iD %*% rnorm(7)
slope <- iD %*% rnorm(7)
#
iD <- t(chol(Vphy))
intercept <- iD %*% rnorm(nspp)
slope <- iD %*% rnorm(nspp)
#
prob <- rep(intercept, each=nsite)
prob <- prob + rep(slope, each=nsite) * rep(env, nspp)
prob <- prob + rnorm(nspp*nsite)
pres <- rbinom(length(prob), size=1, prob=exp(prob)/(1+exp(prob)))
#
site <- factor(rep(1:nsite, nspp))
species <- factor(rep(1:nspp, each=nsite))
env <- rep(env, nspp)
#
r.intercept.spp.indep <- list(1, sp = species, covar = diag(nspp))
r.intercept.spp.phy <- list(1, sp = species, covar = Vphy)
r.slope.spp.indep <- list(env, sp = species, covar = diag(nspp))
r.slope.spp.phy <- list(env, sp = species, covar = Vphy)
r.site <- list(1, site = site, covar = diag(nsite))
rnd.effects <- list(r.intercept.spp.indep, r.intercept.spp.phy, r.slope.spp.indep, r.slope.spp.phy, r.site)
#
model <- communityPGLMM(pres ~ env, family = "binomial", sp = species, site = site, random.effects = rnd.effects, REML = TRUE, verbose = FALSE)
communityPGLMM.binary.LRT(model, re.number = 1)
## $LR
## [1] -1.933817e-05
##
## $df
## [1] 1
##
