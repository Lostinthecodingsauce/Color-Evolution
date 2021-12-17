######Ancestral State recosntruction for life history traits 
#### Reduced Patches Script 
library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("RPANDA")
library(ggplot2)
library(phylosignal)
library(factoextra)
library("evomap")
library(caper)
library(dplyr)
library(jtools)
# https://cran.r-project.org/web/packages/jtools/vignettes/summ.html#Report_robust_standard_errors
library(kableExtra)
library(huxtable)

setwd("C:/Users/bfsco/Desktop/Masters Research/Chapter 2/analysis/WPTSC/MT_tree_analysis")

# Read in Tree and data 
Card.tree <- read.nexus("FullMCC.tree.nexus")
plot(Card.tree)

to.drop <- ("Cyanocompsa_cyanoides")
Card.tree2 <- drop.tip(Card.tree,to.drop)

data <- read.csv("Chap2data.csv", header = TRUE) 


rownames(data) <- data[,1,]
log.Brill <- log(data[c(13,14)])
data<- cbind(data,log.Brill)
data<- data[,c(-13,-14)]
df<- data[,c(-1)]


df<-df[Card.tree$tip.label, ] 


comp.data<-comparative.data(Card.tree, data, names.col = "phylo", vcv=TRUE, warn.dropped=TRUE)
comp.data2<-comparative.data(Card.tree2, data, names.col = "phylo", vcv=TRUE, warn.dropped=TRUE)
df2 <- comp.data2$data
###################
a<-pgls(fPC2~mPC2, data=comp.data2,lambda = "ML") #
summary(a)

with(df2, plot(fPC2~mPC2, xlab = "mPC2", ylab = "fPC2",points(cex = .5, col = "dark red"),
              main = "Female PC2 vs \ Male PC2"))

ggplot(full,aes(y=fPC1,x=Climate.PC2))+geom_point(colour = "blue", size = 5)+ theme_light(base_size=40)

abline(a)
#########################################
a <- ggplot(df2,aes(y=fPC1,x=mPC1))+geom_point(colour = "black",size = 5)+stat_smooth(method="lm",se=FALSE,colour = "black", size = 3)
                 a    + ggtitle("Female PC1 vs Male PC1") + theme_light(base_size=20) 

ggplot(df2,aes(y=fPC2,x=mPC2))+geom_point(colour = "black",size = 5)+stat_smooth(method="lm",se=FALSE,colour = "black", size = 3)+ theme_light(base_size=20)

abline(a)
#########################################
b<-pgls(fAvgBrill~Habitat+Migration, data=comp.data,lambda = "ML") #
summary(b)
plot(b)
abline(b)
b$aicc

with(df2, plot(fPC2~mPC2, xlab = "fPC2", ylab = "mPC2",
                main = "Female PC2 vs \ Male PC2"))
abline(a)







df$fAvgBrill
b<-pgls(fAvgBrill~Migration + Strata, data=comp.data,lambda = "ML") #
summary(b)
b$aicc

plot(b)
with(df2, plot(fPC2~mPC2, xlab = "fPC2", ylab = "mPC2", 
              main = "Female PC2 vs \ Male PC2")) 
abline(b)

library(lme4)
fm1 <- glm(fPC1~mPC1, data = df2, family = gaussian)
summ(fm1)
effect_plot(fm1, pred = mPC1, interval = TRUE, plot.points = TRUE)

b<-pgls(mPC2~Forest.Dependency+Migration, data=comp.data,lambda = "ML") #
summary(b)
b$aicc
plot(b)
abline(b)
b$aicc




######### Plots and t-tests ############
df$mAvgBrill


df %>%
  ggplot( aes(x=nonforest, y=mPC2, fill=nonforest)) + 
  geom_boxplot() +
  xlab("Migration") +
  theme(legend.position="none") +
  xlab("") +theme_light(base_size=20)

sp<-ggplot(data, aes(x=fPC1, y=fPC2)) +  geom_point(size = 4)
sp+ theme_light(base_size=20)


sp3<-ggplot(df, aes(x=mPC1, y=mPC2, color=Forest.Dependency)) + geom_text(label=rownames(df))
sp3+scale_color_gradient(low="blue", high="red")

sp2<-ggplot(df, aes(x=mPC1, y=mPC2)) +  geom_point(aes(colour = HabitatA), size = 3)
sp2+theme_bw() 

# geom_boxplot proposes several arguments to custom appearance
ggplot(df, aes(x=fPC1, y=Strata)) + 
  geom_boxplot(
    # custom boxes
    color="blue",
    fill="blue",
    alpha=0.2,
    # Notch?
    notch=TRUE,
    notchwidth = 0.8) #,
     # custom outliers
   ## outlier.colour="red",
    #outlier.fill="red",
    #outlier.size=3  )


test <- lm(fPC2~mPC2, data = data)
summ(test)

m1FD <- gls(mPC1~Forest.Dependency, data = data)
m2FD <- gls(mPC2~Forest.Dependency, data = data)
f1FD <- gls(fPC1~Forest.Dependency, data = data)
f2FD <- gls(fPC2~Forest.Dependency, data = data)

plot_summs(m1FD, m2FD, f1FD, f2FD,scale = TRUE,plot.distributions = FALSE,
           model.names = c("mPC1", "mPC2", "fPC1","fPC2"))

m1HA <- gls(mPC1~HabitatA, data = data)
m2HA <- gls(mPC2~HabitatA, data = data)
f1HA <- gls(fPC1~HabitatA, data = data)
f2HA <- gls(fPC2~HabitatA, data = data)

plot_summs(m1HA, m2HA, f1HA, f2HA,scale = TRUE,plot.distributions = FALSE,
           model.names = c("mPC1", "mPC2", "fPC1","fPC2"),y.label = "Habitat Metric")

m2FD <- lm(mPC2~Forest.Dependency, data = data)
m2HA <- lm(mPC2~HabitatA, data = data)
m2h <- lm(mPC2~Habitat, data = data)



#####################################################################################3
#Plot ~ Including only significant ones, and all NDVI+ all Tree Cover 
effect_plot(m2FD, pred = Forest.Dependency, interval = TRUE, plot.points = TRUE)


plot_summs(m1_Land1.a, m2_Land1.a, f1_Land1.a, f2_Land1.a,scale = TRUE,plot.distributions = TRUE,
           model.names = c("mPC1", "mPC2", "fPC1","fPC2"))

export_summs(m1_Land1.a, m2_Land1.a, f1_Land1.a, f2_Land1.a, scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]")






chisq.test(df$HabitatA, df$Forest.Dependency, correct=FALSE)

#####################################################################

####################################################################3
## Categorical reconstruction below. First reconstruct cont. traits
Strata<-as.factor(setNames(data[,9],rownames(data)))
StrataA<-as.factor(setNames(data[,10],rownames(data)))
Habitat<-as.factor(setNames(data[,7],rownames(data)))
HabitatA<-as.factor(setNames(data[,8],rownames(data)))
Forest.dep<-as.factor(setNames(data[,11],rownames(data)))
Mig<-as.factor(setNames(data[,12],rownames(data)))
Iris<-as.factor(setNames(data[,15],rownames(data)))
Bill<-as.factor(setNames(data[,16],rownames(data)))

## Dot plots of traits on tree
dotTree(Card.tree,Habitat)
dotTree(Card.tree,HabitatA)
dotTree(Card.tree,Strata)
dotTree(Card.tree,Mig)
dotTree(Card.tree,Forest.dep)
dotTree(Card.tree,Iris)
dotTree(Card.tree,Bill)

Card.tree <- read.nexus("MCC.Cardinalidae.nexus")
### Habitat 
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
cols<-setNames(c("red","blue","green"),levels(HabitatA))
tiplabels(pie=to.matrix(HabitatA[Card.tree$tip.label],
                        levels(HabitatA)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)
#### Ancestral State Recon
fitER<-ace(HabitatA,Card.tree,model="ER",type="discrete")
fitER
fitER$lik.anc #marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities
#Overlay these posterior probabilities on the tree:
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:Card.tree$Nnode+Ntip(Card.tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(Habitat[Card.tree$tip.label],
                        levels(Habitat)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

### Nonforest #####################################################
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
cols<-setNames(c("red","blue"),levels(Habitat))
tiplabels(pie=to.matrix(Habitat[Card.tree$tip.label],
                        levels(Habitat)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

#### Ancestral State Recon
fitER<-ace(Habitat,Card.tree,model="ER",type="discrete")
fitER
fitER$lik.anc #marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities
#Overlay these posterior probabilities on the tree:
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:Card.tree$Nnode+Ntip(Card.tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(Habitat[Card.tree$tip.label],
                        levels(Habitat)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)


### Strata #####################################################
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
cols<-setNames(c("red","blue","green","yellow"),levels(Strata))
tiplabels(pie=to.matrix(Strata[Card.tree$tip.label],
                        levels(Strata)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)
#### Ancestral State Recon
fitER<-ace(Strata,Card.tree,model="ER",type="discrete")
fitER
fitER$lik.anc #marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities
#Overlay these posterior probabilities on the tree:
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:Card.tree$Nnode+Ntip(Card.tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(Strata[Card.tree$tip.label],
                        levels(Strata)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)


### Migration #####################################################
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
cols<-setNames(c("red","blue"),levels(Mig))
tiplabels(pie=to.matrix(Mig[Card.tree$tip.label],
                        levels(Mig)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

#### Ancestral State Recon
fitER<-ace(Mig,Card.tree,model="ER",type="discrete")
fitER
fitER$lik.anc #marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities
#Overlay these posterior probabilities on the tree:
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:Card.tree$Nnode+Ntip(Card.tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(Mig[Card.tree$tip.label],
                        levels(Mig)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

### Forest Dependency
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
cols<-setNames(c("red","blue","green","yellow"),levels(Forest.dep))
tiplabels(pie=to.matrix(Forest.dep[Card.tree$tip.label],
                        levels(Forest.dep)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)
#### Ancestral State Recon
fitER<-ace(Forest.dep,Card.tree,model="ER",type="discrete")
fitER
fitER$lik.anc #marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities
#Overlay these posterior probabilities on the tree:
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:Card.tree$Nnode+Ntip(Card.tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(Forest.dep[Card.tree$tip.label],
                        levels(Forest.dep)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

###################
### Bill #####################################################
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
cols<-setNames(c("red","blue"),levels(Bill))
tiplabels(pie=to.matrix(Bill[Card.tree$tip.label],
                        levels(Bill)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

#### Ancestral State Recon
fitER<-ace(Bill,Card.tree,model="ER",type="discrete")
fitER
fitER$lik.anc #marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities
#Overlay these posterior probabilities on the tree:
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:Card.tree$Nnode+Ntip(Card.tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(Bill[Card.tree$tip.label],
                        levels(Bill)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

### Iris #####################################################
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
cols<-setNames(c("red","blue"),levels(Iris))
tiplabels(pie=to.matrix(Iris[Card.tree$tip.label],
                        levels(Iris)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

#### Ancestral State Recon
fitER<-ace(Iris,Card.tree,model="ER",type="discrete")
fitER
fitER$lik.anc #marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities
#Overlay these posterior probabilities on the tree:
plotTree(Card.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:Card.tree$Nnode+Ntip(Card.tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(Iris[Card.tree$tip.label],
                        levels(Iris)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)



#####

### MCMC reconstruction
mtree<-make.simmap(Card.tree,Mig,model="ER")
plot(mtree,cols,type="fan",fsize=0.7,ftype="i")
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

mtrees<-make.simmap(Card.tree,Mig,model="ER",nsim=100)
par(mfrow=c(10,10))
null<-sapply(mtrees,plotSimmap,colors=cols,lwd=1,ftype="off")

pd<-summary(mtrees)
pd

plot(pd,fsize=0.6,ftype="i",colors=cols,ylim=c(-2,Ntip(Card.tree)))
add.simmap.legend(colors=cols[2:1],prompt=FALSE,x=0,y=-4,vertical=FALSE)

plot(sample(mtrees,1)[[1]],cols,fsize=0.6,ftype="i",
     ylim=c(-2,Ntip(Card.tree)))
nodelabels(pie=pd$ace,piecol=cols,cex=0.5)
add.simmap.legend(colors=cols[2:1],prompt=FALSE,x=0,y=-4,
                  vertical=FALSE)