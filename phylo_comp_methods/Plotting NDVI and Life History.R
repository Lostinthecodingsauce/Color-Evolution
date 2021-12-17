### Plotting Data 
##################################################################

setwd("C:/Users/bfsco/Desktop/Masters Research/Plumage/analysis/WPTSC/MT_tree_analysis")
#tip <- c("Cyanocompsa_cyanoides") #"Cyanocompsa_cyanoides"
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- read.nexus("FullMCC.tree.nexus")
#Card.tree <- drop.tip(FullCard.tree,tip)
library("ggpubr")
library(phytools)

#All <- read.csv("MCMC_NDVI.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') ##Raw Spatial Data
newer <- read.csv("new.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') # PCA not log-transformed 

raw.traits <- read.csv("RawTraits_NDVI.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') # Non-scaled spatial data
raw.traits <- subset(raw.traits[c(2:21)]) # get only raw traits 
# Transform Spatial data to have correct scale 
df.mean <- All[c(16,21,26,31)]*0.0001
raw.traits <- subset(raw.traits[c(2:21)])




# Transform Spatial data to have correct scale based on NDVI values 
df.MODIS <- newer[c(16:21)]*0.0001 # df.mean <- newer[c(16,21,26,31)]*0.0001
df.GPP <- newer[c(22:24)]*0.01 # df.mean <- newer[c(16,21,26,31)]*0.0001

temp <- subset(newer[c(1:15)]) # subset(All[c(1:15,17:20,22:25,27:30,32:36)])
df <- cbind.data.frame(temp,df.MODIS)
data <- cbind.data.frame(df,df.GPP)

#####

# If need species to be row.names (for phylo purposes)

####Make species header for rows
full <-data
rownames(full) <- full[,1,]
full<- full[-c(1)]

####Make species header for rows
rownames(raw.traits) <- raw.traits[,1,]
raw.traits<- raw.traits[-c(1)]


################################################################

## Create martix for each PC score
fuck <- subset.data.frame(full[c(1:7,14:23)])

Tree.Cover<-as.matrix(fuck)[,7]
NDVI.mean<-as.matrix(fuck)[,9]
NIR.mean<-as.matrix(fuck)[,13]
NDVI.stdDev<-as.matrix(fuck)[,10]
NIR.stdDev<-as.matrix(fuck)[,14]
NIR.Var<-as.matrix(fuck)[,17]

##### Plotting ################
lims<-c(floor(1),ceiling(100))
obj.Tree.Cover<-contMap(Card.tree,Tree.Cover, fsize=c(0.6,1),outline=FALSE, legend= TRUE)

### Invert color of tree
Tree<-setMap(obj.Tree.Cover, invert = TRUE)
plot(Tree)


### mess around with colors 
pal <- colorRampPalette(c("red", "green"))
pal(20)

Tree<-setMap(obj.Tree.Cover, pal(10))
plot(Tree)



obj.NIR.mean<-contMap(Card.tree,NIR.mean, fsize=c(0.6,1),outline=FALSE)
obj.NDVI.mean<-contMap(Card.tree,NDVI.mean, fsize=c(0.6,1),outline=FALSE)
obj.NIR.stdDev<-contMap(Card.tree,NIR.stdDev, fsize=c(0.6,1),outline=FALSE)
obj.NIR.Var<-contMap(Card.tree,NIR.Var, fsize=c(0.6,1),outline=FALSE)


obj.NDVI.stdDev<-contMap(Card.tree,NDVI.stdDev, fsize=c(0.6,1),outline=FALSE)



##########################################################################################

fmode<-as.factor(setNames(full[,8],rownames(full)))

dotTree(Card.tree,fmode,colors=setNames(c("blue","red"),
                                        c("open","closed")),ftype="i",fsize=0.7)
text(x=28,y=-1.8,"Habitat",pos=4)
plot(obj.Tree.Cover$tree,colors=obj.Tree.Cover$cols,add=TRUE,lwd=6,ftype="off",outline=TRUE,fsize=c(0,1.2),
    xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)
#add.color.bar(leg = c(0,100), cols=obj.Tree.Cover$cols,title="Percent.Tree.cover",outline = FALSE),
 #             fsize=c(0,1.2),prompt=TRUE,x=0,y=100)

obj<-contMap(Card.tree,Tree.Cover, fsize=c(0.6,1),outline=FALSE, legend= TRUE)

h<-max(Percent.Tree.cover(obj$tree))
plot(obj,legend=FALSE,xlim=c(-0.25*h,1.1*h))
## now add our legend, but without the text underneath
add.color.bar(Ntip(obj$tree)-3,obj$cols,title="trait value",
              lims=NULL,digits=3,direction="upwards",
              subtitle="",lwd=15,x=-0.2,y=2,prompt=FALSE)
## get line width in user units
LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
lines(x=rep(-0.2*h+LWD*15/2,2),y=c(2,Ntip(obj$tree)-1))
nticks<-10
Y<-cbind(seq(2,Ntip(obj$tree)-1,length.out=nticks),
         seq(2,Ntip(obj$tree)-1,length.out=nticks))
X<-cbind(rep(-0.2*h+LWD*15/2,nticks),
         rep(-0.2*h+LWD*15/2+0.02*h,nticks))
for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
ticks<-seq(obj$lims[1],obj$lims[2],length.out=nticks)
text(x=X[,2],y=Y[,2],round(ticks,3),pos=4,cex=0.8)


#################################
fmode<-as.factor(setNames(full[,9],rownames(full)))

dotTree(Card.tree,fmode,colors=setNames(c("blue","red","green"),
                                        c("O","C","U")),ftype="i",fsize=0.7)
text(x=28,y=-1.8,"Strata",pos=4)
plot(obj.Tree.Cover$tree,colors=obj.Tree.Cover$cols,add=TRUE,lwd=6,ftype="off",outline=TRUE,fsize=c(0,1.2),
     xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)
############################################################################################################
############################################################################################################

############################################################################################################
############################################################################################################



#Plot ~ Including only significant ones, and all NDVI+ all Tree Cover 
effect_plot(NIR.mean_m1, pred = NIR.mean, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.mean_m2, pred = NIR.mean, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.mean_f1, pred = NIR.mean, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.mean_f2, pred = NIR.mean, interval = TRUE, plot.points = TRUE)

## VisRed
effect_plot(VisRed.mean_m1, pred = VisRed.mean, interval = TRUE, plot.points = TRUE)
effect_plot(VisRed.mean_m2, pred = VisRed.mean, interval = TRUE, plot.points = TRUE)
effect_plot(VisRed.mean_f1, pred = VisRed.mean, interval = TRUE, plot.points = TRUE)
effect_plot(VisRed.mean_f2, pred = VisRed.mean, interval = TRUE, plot.points = TRUE)

## Tree Cover 
effect_plot(f1tree, pred = Percent.Tree.cover, interval = TRUE, plot.points = TRUE)
effect_plot(m1tree, pred = Percent.Tree.cover, interval = TRUE, plot.points = TRUE)
effect_plot(f2tree, pred = Percent.Tree.cover, interval = TRUE, plot.points = TRUE)
effect_plot(m2tree, pred = Percent.Tree.cover, interval = TRUE, plot.points = TRUE)


## NDVI.mean
effect_plot(NDVI.mean_m1, pred = NDVI.mean, interval = TRUE, plot.points = TRUE)
effect_plot(NDVI.mean_m2, pred = NDVI.mean, interval = TRUE, plot.points = TRUE)
effect_plot(NDVI.mean_f1, pred = NDVI.mean, interval = TRUE, plot.points = TRUE)
effect_plot(NDVI.mean_f2, pred = NDVI.mean, interval = TRUE, plot.points = TRUE)

# NIR.stdDev
## NIR.stdDev
effect_plot(NIR.stdDev_m1, pred = NIR.stdDev, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.stdDev_m2, pred = NIR.stdDev, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.stdDev_f1, pred = NIR.stdDev, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.stdDev_f2, pred = NIR.stdDev, interval = TRUE, plot.points = TRUE)

## NDVI.stdDev
effect_plot(NDVI.stdDev_m1, pred = NDVI.stdDev, interval = TRUE, plot.points = TRUE)
effect_plot(NDVI.stdDev_m2, pred = NDVI.stdDev, interval = TRUE, plot.points = TRUE)
effect_plot(NDVI.stdDev_f1, pred = NDVI.stdDev, interval = TRUE, plot.points = TRUE)
effect_plot(NDVI.stdDev_f2, pred = NDVI.stdDev, interval = TRUE, plot.points = TRUE)

## NIR.Var
effect_plot(NIR.Var_m1, pred = NIR.Var, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.Var_m2, pred = NIR.Var, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.Var_f1, pred = NIR.Var, interval = TRUE, plot.points = TRUE)
effect_plot(NIR.Var_f2, pred = NIR.Var, interval = TRUE, plot.points = TRUE)

# Fuck yes great plot S
plot_summs(NIR.mean_m1, NIR.mean_m2,NIR.mean_f1, NIR.mean_f2,scale = TRUE,plot.distributions = TRUE,
           model.names = c("mPC1", "mPC2", "fPC1","fPC2"))

export_summs(NIR.mean_m1, NIR.mean_m2,NIR.mean_f1, NIR.mean_f2, scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]")

plot_summs(NDVI.mean_m1, NDVI.mean_m2,NDVI.mean_f1, NDVI.mean_f2,scale = TRUE,plot.distributions = TRUE,
           model.names = c("mPC1", "mPC2", "fPC1","fPC2"))

export_summs(NDVI.mean_m1, NDVI.mean_m2,NDVI.mean_f1, NDVI.mean_f2, scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]")

plot_summs(m1tree, m2tree,f1tree, f2tree,scale = TRUE,plot.distributions = TRUE,
           model.names = c("mPC1", "mPC2", "fPC1","fPC2"))


###################################################################################









############################
All <- data

library("ggpubr")

ggscatter(All, x = "NDVI.mean", y = "mPC1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "NDVI.mean~mPC1",
          xlab = "NDVI.mean", ylab = "mPC1")

ggscatter(All, x = "NDVI.mean", y = "mPC2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "NDVI.mean~mPC2",
          xlab = "NDVI.mean", ylab = "mPC2")

ggscatter(All, x = "NDVI.mean", y = "fPC1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "NDVI.mean~fPC1",
          xlab = "NDVI.mean", ylab = "fPC1")

ggscatter(All, x = "NDVI.mean", y = "fPC2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "NDVI.mean~fPC2",
          xlab = "NDVI.mean", ylab = "fPC2")

####
ggscatter(All, x = "Percent.Tree.Cover", y = "mPC1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "Percent.Tree.Cover~mPC1",
          xlab = "Percent.Tree.Cover", ylab = "mPC1")

ggscatter(All, x = "Percent.Tree.Cover", y = "mPC2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "Percent.Tree.Cover~mPC2",
          xlab = "Percent.Tree.Cover", ylab = "mPC2")

ggscatter(All, x = "Percent.Tree.Cover", y = "fPC1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "Percent.Tree.Cover~fPC1",
          xlab = "Percent.Tree.Cover", ylab = "fPC1")

ggscatter(full, x = "Percent.Tree.Cover", y = "VisRed.mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "Percent.Tree.Cover~VisRed.mean",
          xlab = "Percent.Tree.Cover", ylab = "VisRed.mean")


# Scatter plot


sp2<-ggplot(full, aes(x=Percent.Tree.cover, y=mPC1, color=NDVI.stdDev)) + geom_point(size = 3)
#sp3<-ggplot(full, aes(x=mPC1, y=Percent.Tree.cover, color=NDVI.stdDev)) + geom_text(label=rownames(full))
sp2+scale_color_gradient(low="blue", high="red")



sp3<-ggplot(full, aes(x=fPC2, y=Percent.Tree.cover, color=VisRed.stdDev)) + geom_text(label=rownames(full))
sp3+scale_color_gradient(low="blue", high="red")

sp4<-ggplot(full, aes(x=mPC1, y=mPC2, color=VisRed.mean)) + geom_text(label=rownames(full))
sp4+scale_color_gradient(low="blue", high="red")




###
# cLUSTERING 
#### 
cluster <- All[-c(1,8:13)]
head(cluster)
cls <- kmeans(x = cluster, centers = 5)
cluster$cluster <- as.character(cls$cluster)
head(cluster)

ggplot() + geom_point(data = cluster, 
                   mapping = aes(x = Percent.Tree.cover, 
                                 y = mPC2, 
                                 colour = cluster),size = 3)

############### Boxplots ###############
boxplot(Percent.Tree.cover~Forest.Dependency,data=data, main="Percent.Tree.cover vs Forest Dependency",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Forest", ylab="Percent.Tree.cover")
boxplot(Percent.Tree.cover~Strata,data=data, main="Percent.Tree.cover vs Strata",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Strata", ylab="Percent.Tree.cover")
boxplot(Percent.Tree.cover~Habitat,data=data, main="Percent.Tree.cover vs Habitat",
        col=(c("darkgreen","red")),
        xlab="Habitat", ylab="Percent.Tree.cover")
boxplot(Percent.Tree.cover~Migration,data=data, main="Percent.Tree.cover vs Mig",
        col=(c("darkgreen","red")),
        xlab="Migration", ylab="Percent.Tree.cover")


###

boxplot(NDVI.mean~Forest.Dependency,data=data, main="meanNDVI vs Forest Dependency",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Forest", ylab="mean.NDVI")
boxplot(NDVI.mean~Strata,data=data, main="meanNDVI vs Strata",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Strata", ylab="mean.NDVI")
boxplot(NDVI.mean~Habitat,data=data, main="meanNDVI vs Habitat",
        col=(c("darkgreen","red")),
        xlab="Habitat", ylab="mean.NDVI")
boxplot(NDVI.mean~Migration,data=data, main="meanNDVI vs Mig",
        col=(c("darkgreen","red")),
        xlab="Migration", ylab="mean.NDVI")

boxplot(NDVI.stdDev~Forest.Dependency,data=data, main="stdDevNDVI vs Forest Dependency",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Forest", ylab="stdDev.NDVI")
boxplot(NDVI.stdDev~Strata,data=data, main="stdDevNDVI vs Strata",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Strata", ylab="stdDev.NDVI")
boxplot(NDVI.stdDev~Habitat,data=data, main="stdDevNDVI vs Habitat",
        col=(c("darkgreen","red")),
        xlab="Habitat", ylab="stdDev.NDVI")
boxplot(NDVI.stdDev~Migration,data=data, main="stdDevNDVI vs Mig",
        col=(c("darkgreen","red")),
        xlab="Migration", ylab="stdDev.NDVI")
########################


#######################################################################
### NIR ####
boxplot(NIR.mean~Forest.Dependency,data=data, main="meanNIR vs Forest Dependency",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Forest", ylab="mean.NIR")
boxplot(NIR.mean~Strata,data=data, main="meanNIR vs Strata",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Strata", ylab="mean.NIR")
boxplot(NIR.mean~Habitat,data=data, main="meanNIR vs Habitat",
        col=(c("darkgreen","red")),
        xlab="Habitat", ylab="mean.NIR")
boxplot(NIR.mean~Migration,data=data, main="meanNIR vs Mig",
        xlab="Migration", ylab="mean.NIR")

boxplot(NIR.stdDev~Forest.Dependency,data=data, main="stdDevNIR vs Forest Dependency",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Forest", ylab="stdDev.NIR")
boxplot(NIR.stdDev~Strata,data=data, main="stdDevNIR vs Strata",
        col=(c("gold","darkgreen","blue","red")),
        xlab="Strata", ylab="stdDev.NIR")
boxplot(NIR.stdDev~Habitat,data=data, main="stdDevNIR vs Habitat",
        col=(c("darkgreen","red")),
        xlab="Habitat", ylab="stdDev.NIR")
boxplot(NIR.stdDev~Migration,data=data, main="stdDevNIR vs Mig",
        xlab="Migration", ylab="stdDev.NIR")

###################################################



hist(data$NIR.mean)
hist(data$NIR.stdDev)

hist(data$NDVI.mean)
hist(data$NDVI.stdDev)


hist(data$GPP.mean)
hist(data$GPP.stdDev)

hist(data$Percent.Tree.cover)

###############

library(ggimage)
library(ggtree)
library(ggplot2)


set.seed(2019-10-31)
tr <- rtree(10)

d1 <- data.frame(
        # only some labels match
        label = c(FullCard.tree$tip.label[sample(5, 5)], "NDVI"),
        value = sample(1:6, 6))

d2 <- data.frame(
        label = rep(FullCard.tree$tip.label,10),
        category = rep(LETTERS[1], each=1),
        value = rnorm(44, 0, 3)) 

g <- ggtree(FullCard.tree) + geom_tiplab(align=TRUE)

p1 <- ggplot(d1, aes(label, value)) + geom_col(aes(fill=label)) + 
        geom_text(aes(label=label, y= value+.1)) +
        coord_flip() + theme_tree2() + theme(legend.position='none')

p2 <- ggplot(d2, aes(x=category, y=label)) + 
        geom_tile(aes(fill=value)) + scale_fill_viridis_c() + 
        theme_tree2() 

cowplot::plot_grid(g, p2, p1, ncol=3) 

library(aplot)
p2 %>% insert_left(g) %>% insert_right(p1, width=.5) 

##################################
remote_folder <- paste0("https://raw.githubusercontent.com/katholt/",
                         "plotTree/master/tree_example_april2015/")
#remote_folder <- "data/tree_example_april2015/" 

## read the phylogenetic tree
tree <- read.tree(paste0(remote_folder, "tree.nwk"))

## read the sampling information data set
info <- read.csv(paste0(remote_folder,"info.csv"))

## read and process the allele table
snps<-read.csv(paste0(remote_folder, "alleles.csv"), header = F,
               row.names = 1, stringsAsFactor = F)
snps_strainCols <- snps[1,] 
snps<-snps[-1,] # drop strain names
colnames(snps) <- snps_strainCols

gapChar <- "?"
snp <- t(snps)
lsnp <- apply(snp, 1, function(x) {
        x != snp[1,] & x != gapChar & snp[1,] != gapChar
})
lsnp <- as.data.frame(lsnp)
lsnp$pos <- as.numeric(rownames(lsnp))
lsnp <- tidyr::gather(lsnp, name, value, -pos)
snp_data <- lsnp[lsnp$value, c("name", "pos")]

## read the trait data
bar_data <- read.csv(paste0(remote_folder, "bar.csv"))

## visualize the tree 
p <- ggtree(tree) 
q <- ggtree(FullCard.tree) 
## attach the sampling information data set 
## and add symbols colored by location
p <- p %<+% info + geom_tippoint(aes(color=location))
q1 <- q %<+% full + geom_tippoint(aes(color=Habitat))
## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p + geom_facet(panel = "SNP", data = snp_data, geom = geom_point, 
               mapping=aes(x = pos, color = location), shape = '|') +
        geom_facet(panel = "Trait", data = bar_data, geom = ggstance::geom_barh, 
                   aes(x = dummy_bar_value, color = location, fill = location), 
                   stat = "identity", width = .6) +
        theme_tree2(legend.position=c(.05, .85))






p <- ggtree(Card.tree) 
p1 <- p %<+% data + geom_tippoint(aes(color=Habitat))
##################################################################




## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p1 + geom_facet(panel = "Percent.Tree.Cover", data = Forest.Cover, geom = ggstance::geom_barh, 
                   aes(x = Percent.Tree.Cover, color = green), 
                   stat = "identity", width = .6) +
        theme_tree2(legend.position=c(.05, .85))


###################################################
### Violin PLot ############

# Violin Plots
library(vioplot)

x1 <- All$NDVI.mean[All$Forest.Dependency=="High"]
x2 <- All$NDVI.mean[All$Forest.Dependency=="Medium"]
x3 <- All$NDVI.mean[All$Forest.Dependency=="Low"]
x4 <- All$NDVI.mean[All$Forest.Dependency=="None"]
vioplot(x1, x2, x3, x4, names=c("High", "Medium", "Low","None"),
        col="gold")
title("Violin Plots of Mean.NDVI~Forest.Dependency")



## Chi=squared 
library('MASS')
testa1 <- chisq.test(All$NDVI.mean, All$Forest.Dependency)
testa2 <- chisq.test(All$NDVI.stdDev, All$Forest.Dependency)
testb1 <- chisq.test(All$NDVI.mean, All$StrataA)
testb2 <- chisq.test(All$NDVI.stdDev, All$StrataA)
testc1 <- chisq.test(All$NDVI.mean, All$Mig)
testc2 <- chisq.test(All$NDVI.stdDev, All$Mig)
testd1 <- chisq.test(All$NDVI.mean, All$Habitat)
testd2 <- chisq.test(All$NDVI.stdDev, All$Habitat)

### Correlated Chi-squared
library('DescTools')


cora1 <- ContCoef(All$NDVI.mean, All$Forest.Dependency)
cora2 <- ContCoef(All$NDVI.stdDev, All$Forest.Dependency)
corb1 <- ContCoef(All$NDVI.mean, All$StrataA)
corb2 <- ContCoef(All$NDVI.stdDev, All$StrataA)
corc1 <- ContCoef(All$NDVI.mean, All$Mig)
corc2 <- ContCoef(All$NDVI.stdDev, All$Mig)
cord1 <- ContCoef(All$NDVI.mean, All$Habitat)
cord2 <- ContCoef(All$NDVI.stdDev, All$Habitat)
