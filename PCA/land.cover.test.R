############ LAND COVER PCA TEST #######################



library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("RPANDA")
library(ggplot2)
library(phylosignal)
library(factoextra)
library("ggpubr")
library(jtools)
# https://cran.r-project.org/web/packages/jtools/vignettes/summ.html#Report_robust_standard_errors
library(kableExtra)
library(huxtable)

setwd("C:/Users/bfsco/Desktop/Masters Research/Plumage/analysis/WPTSC/MT_tree_analysis")


 # read in data 
data <- read.csv("landCover.test.csv", header = TRUE) # values are divided by range, = prortion of each land cover in range
#data <- read.csv("Climate.SD.csv", header = TRUE) 

####Make species header for rows
rownames(data) <- data[,1,]
data<- data[c(-1)]



df <- round(data, digits = 2)

test<- df[-c(1,3,6,11,13,15)]


#### 1) Regular PCA   ##############################################################################

pca <- princomp(df, cor = FALSE, scores = TRUE)
pca <- princomp(df,cor = TRUE, scores = TRUE)

summary(pca)#Importance of Components
print(pca) #Standard Deviations

## Get scores
scores <- as.data.frame(pca$scores)
pca$loadings

fviz_eig(pca)



fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

write.csv(scores, file = "landcoverPCAscores.csv")

####################### GLS ########################################

data <- read.csv("new.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') 

#data <- read.csv("test.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') 


### Change to appropriate scale for spatial data. 
df.toscaled <- data[c(16:21)]*0.0001

temp <- subset(data[c(1:15,22:32)])
data <- cbind.data.frame(temp,df.toscaled) #df.toscale
###########################################################
####Make species header for rows
full <-data
rownames(full) <- full[,1,]
full<- full[-c(1)]


##################################################################################
#### Plotting ####
##################################################################################
library(jtools)
# https://cran.r-project.org/web/packages/jtools/vignettes/summ.html#Report_robust_standard_errors
library(kableExtra)
library(huxtable)



# .1
m1_LandComp.1<-lm(mPC1~LandComp.1, data=data) 
m2_LandComp.1<-lm(LandComp.1~mPC2, data=data) # # significant 
f1_LandComp.1<-lm(fPC1~LandComp.1, data=data) # significant 
f2_LandComp.1<-lm(fPC2~LandComp.1, data=data)

summ(m1_LandComp.1,  digits = 5)
summ(m2_LandComp.1,  digits = 5)
summ(f1_LandComp.1,  digits = 5)
summ(f2_LandComp.1,  digits = 5)



# Comp.2
m1_LandComp.2<-lm(mPC1~LandComp.2, data=data) 
m2_LandComp.2<-lm(mPC2~LandComp.2, data=data)
f1_LandComp.2<-lm(fPC1~LandComp.2, data=data) 
f2_LandComp.2<-lm(fPC2~LandComp.2, data=data)

summ(m1_LandComp.2,  digits = 5)
summ(m2_LandComp.2,  digits = 5)
summ(f1_LandComp.2,  digits = 5)
summ(f2_LandComp.2,  digits = 5)



# Comp.3
m1_LandComp.3<-lm(mPC1~LandComp.3, data=data) 
m2_LandComp.3<-lm(mPC2~LandComp.3, data=data)
f1_LandComp.3<-lm(fPC1~LandComp.3, data=data) 
f2_LandComp.3<-lm(fPC2~LandComp.3, data=data)

summ(m1_LandComp.3,  digits = 5)
summ(m2_LandComp.3,  digits = 5)
summ(f1_LandComp.3,  digits = 5)
summ(f2_LandComp.3,  digits = 5)

sp<-ggplot(full, aes(x=fPC1, y=fPC2, color=LandComp.1)) + geom_text(label=rownames(full))
sp+scale_color_gradient(low="red", high="darkgreen")

sp2<-ggplot(full, aes(x=LandComp.1, y=mPC2, color=Percent.Tree.cover)) + geom_text(label=rownames(full))
sp2+scale_color_gradient(low="red", high="darkgreen")

sp3<-ggplot(full, aes(x=LandComp.1, y=fPC1, color=Percent.Tree.cover)) + geom_text(label=rownames(full))
sp3+scale_color_gradient(low="red", high="darkgreen")

sp4<-ggplot(full, aes(x=LandComp.1, y=LandComp.2, color=Percent.Tree.cover)) + geom_text(label=rownames(full))
sp4+scale_color_gradient(low="red", high="darkgreen")


sp4<-ggplot(full, aes(x=LandComp.1, y=LandComp.2)) + geom_text(label=rownames(full))
sp4
+scale_color_gradient(low="red", high="darkgreen")


sp3<-ggplot(full, aes(x=Percent.Tree.cover, y=fPC1, color=Habitat)) + geom_text(label=rownames(full))
sp3 
####
pal <- colorRampPalette(c("red", "green"))

sp<-ggplot(full, aes(x=mPC2, y=NIR.mean, color=LandComp.1)) + geom_text(label=rownames(full))
sp+scale_color_gradient(low="darkgreen", high="red")

pal <- colorRampPalette(c("red", "green"))
sp3+scale_colour_brewer(palette = "RdYlBu")

sp3+scale_colour_brewer(palette = "RdYlBu")
scale_fill_brewer
######################################################


#Plot ~ Including only significant ones, and all NDVI+ all Tree Cover 
effect_plot(m1_LandComp.1, pred = LandComp.1, interval = TRUE, plot.points = TRUE)
effect_plot(m2_LandComp.1, pred = LandComp.1, interval = TRUE, plot.points = TRUE)
effect_plot(f1_LandComp.1, pred = LandComp.1, interval = TRUE, plot.points = TRUE)
effect_plot(f2_LandComp.1, pred = LandComp.1, interval = TRUE, plot.points = TRUE)


effect_plot(m1_LandComp.2, pred = LandComp.2, interval = TRUE, plot.points = TRUE)
effect_plot(m2_LandComp.2, pred = LandComp.2, interval = TRUE, plot.points = TRUE)
effect_plot(f1_LandComp.2, pred = LandComp.2, interval = TRUE, plot.points = TRUE)
effect_plot(f2_LandComp.2, pred = LandComp.2, interval = TRUE, plot.points = TRUE)


effect_plot(m1_LandComp.3, pred = LandComp.3, interval = TRUE, plot.points = TRUE)
effect_plot(m2_LandComp.3, pred = LandComp.3, interval = TRUE, plot.points = TRUE)
effect_plot(f1_LandComp.3, pred = LandComp.3, interval = TRUE, plot.points = TRUE)
effect_plot(f2_LandComp.3, pred = LandComp.3, interval = TRUE, plot.points = TRUE)


##########################
### Testing 

m2_Land.1<-glm(mPC2~LandComp.1+NDVI.mean, data=data) # # significant
summ(m2_Land.1,  digits = 5)
effect_plot(m2_Land.1, pred = LandComp.1, interval = TRUE, plot.points = TRUE)






m2_Land.2<-glm(mPC2~NIR.mean+LandComp.1, data=data) # # significant
summ(m2_Land.2,  digits = 5)
effect_plot(m2_Land.2, pred = Forest.Dependency, interval = TRUE, plot.points = TRUE)

sp<-ggplot(full, aes(x=mPC2, y=NDVI.mean, color=LandComp.1)) + geom_text(label=rownames(full))
sp+scale_color_gradient(low="darkgreen", high="red")



m2_Land.4<-glm(mPC2~LandComp.1+Forest.Dependency+Strata, data=data) # # significant
summ(m2_Land.4,  digits = 5)
effect_plot(m2_Land.4, pred = LandComp.1, interval = TRUE, plot.points = TRUE)


plot_summs(m2_Land.1, m2_Land.2,m2_tree,m2_Land.4,scale = TRUE,plot.distributions = TRUE,
           model.names = c("LandPC1", "LandPC1+FD", "m2_tree","LandPC1+FD+Stata"))

export_summs(NIR.mean_m1, NIR.mean_m2,NIR.mean_f1, NIR.mean_f2, scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]")









#####################################
### On Tree
####################################
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- read.nexus("FullMCC.tree.nexus")
library(phytools)


## Create martix for each PC score
fuck <- subset.data.frame(full[c(1:7,15:26)])

Tree.Cover<-as.matrix(fuck)[,7]
Land1<-as.matrix(fuck)[,8]
Land2<-as.matrix(fuck)[,9]
Land3<-as.matrix(fuck)[,10]


NDVI.mean<-as.matrix(fuck)[,11]
NIR.mean<-as.matrix(fuck)[,14]
NDVI.stdDev<-as.matrix(fuck)[,12]
NIR.stdDev<-as.matrix(fuck)[,15]

####################################
pal <- colorRampPalette(c("brown", "green"))


## Plotting

lims<-c(floor(1),ceiling(100))
obj.Tree.Cover<-contMap(Card.tree,Tree.Cover, fsize=c(0.6,1),outline=FALSE, legend= TRUE)

### Invert color of tree
Tree<-setMap(obj.Tree.Cover, invert = TRUE)
plot(Tree)


### mess around with colors 
pal <- colorRampPalette(c("red", "green"))


Tree<-setMap(obj.Tree.Cover, pal(10))
plot(Tree)

library(RColorBrewer)
cols <- brewer.pal(9, "BuGn")
pal <- colorRampPalette(cols)

Tree<-setMap(obj.Tree.Cover, pal(20))
plot(Tree)


##### Other plots 
pal <- colorRampPalette(c("green", "red"))


obj.Land1<-contMap(Card.tree,Land1, fsize=c(0.6,1),outline=FALSE)
L1<-setMap(obj.Land1, pal(10))
plot(L1)


pal <- colorRampPalette(c("red", "green"))
obj.Land2<-contMap(Card.tree,Land2, fsize=c(0.6,1),outline=FALSE)
L2<-setMap(obj.Land2, pal(10))
plot(L2)

obj.Land3<-contMap(Card.tree,Land3, fsize=c(0.6,1),outline=FALSE)
L3<-setMap(obj.Land3, pal(10))
plot(L3)
####################################
### PGLS 
######################################

library(ape)
library(phytools)
library(Rcpp)
library(geiger)
library(RPANDA)
library(caper)
library(MCMCglmm)
library(ggplot2)

comp.data<-comparative.data(Card.tree, data, names.col = "X", vcv=TRUE, warn.dropped=TRUE)

m1tree<-pgls(mPC1~Percent.Tree.cover, data=comp.data,lambda = "ML") #
m2tree<-pgls(mPC2~Percent.Tree.cover, data=comp.data,lambda = "ML")#
f1tree<-pgls(fPC1~Percent.Tree.cover, data=comp.data,lambda = "ML") 
f2tree<-pgls(fPC2~Percent.Tree.cover, data=comp.data,lambda = "ML")

summary(m1tree)
summary(m2tree)
summary(f1tree) # yes 
summary(f2tree)


m1Land1<-pgls(mPC1~LandComp.1, data=comp.data,lambda = "ML") #
m2Land1<-pgls(mPC2~LandComp.1, data=comp.data,lambda = "ML")#
f1Land1<-pgls(fPC1~LandComp.1, data=comp.data,lambda = "ML") 
f2Land1<-pgls(fPC2~LandComp.1, data=comp.data,lambda = "ML")

summary(m1Land1)
summary(m2Land1)
summary(f1Land1) # yes 
summary(f2Land1)

m1Land2<-pgls(mPC1~LandComp.2, data=comp.data,lambda = "ML") #
m2Land2<-pgls(mPC2~LandComp.2, data=comp.data,lambda = "ML")#
f1Land2<-pgls(fPC1~LandComp.2, data=comp.data,lambda = "ML") 
f2Land2<-pgls(fPC2~LandComp.2, data=comp.data,lambda = "ML")

summary(m1Land2)
summary(m2Land2)
summary(f1Land2) # yes 
summary(f2Land2)


#######################
m1Land1<-lm(mPC1~LandComp.1 , data=data) 
m2Land1<-lm(mPC2~LandComp.1, data=data)
f1Land1<-lm(fPC1~LandComp.1 , data=data) # significant 
f2Land1<-lm(fPC2~LandComp.1 , data=data)
summ(m1Land1,  digits = 5)
summ(m2Land1,  digits = 5)
summ(f1Land1,  digits = 5)
summ(f2Land1,  digits = 5)

effect_plot(f1Land1, pred = LandComp.1, interval = TRUE, plot.points = TRUE)

plot_summs(m1Land1, m2Land1,f1Land1, f2Land1,scale = TRUE,plot.distributions = TRUE,
           model.names = c("mPC1Land1", "mPC2Land1", "fPC1Land1","fPC2Land1"))
