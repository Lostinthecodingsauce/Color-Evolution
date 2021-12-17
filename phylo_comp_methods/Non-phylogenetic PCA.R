# Non-phylogenetic PCA

library("phylopath")
library(phytools)
library("ape")
library(Rcpp)
library("geiger")
library("caper")
library(tidyverse)
library(ggplot2)
library("RPANDA")
library(factoextra)
library(RPANDA)
library(dplyr)
library(tidyverse)

setwd("C:/Users/bfsco/Desktop/Masters Research/Chapter 2/analysis/WPTSC/MT_tree_analysis")


WPTCS <- read.csv("WPTCS-MT.csv", header = TRUE)
rownames(WPTCS) <- WPTCS[,1,]
WPTCS<- WPTCS[c(-1)]
log.WPTCS <- log(WPTCS[,1:9])
# Remove "Cyanoloxia_cyanoides")
log.WPTCS.c <- log.WPTCS[-c(21), ] 
# Remove  "Piranga_hepatica"
#log.WPTCS.c <- log.WPTCS[-c(21,22,79,80,85,86), ] 

pca <- princomp(log.WPTCS, cor = TRUE, scores = TRUE)
pca.c <- princomp(log.WPTCS.c, cor = TRUE, scores = TRUE)
summary(pca)#Importance of Components
print(pca) #Standard Deviations

## Get scores
scores <- as.data.frame(pca$scores)
pca$loadings
pca.c$loadings



All.scores <- scores %>% 
  rename(
    AllPC.1 = Comp.1,
    AllPC.2 = Comp.2,
    AllPC.3 = Comp.3)


All.scores <- All.scores[c(1:3)]
All.scores$species <- row.names(All.scores)
sp<-ggplot(All.scores, aes(x=AllPC.1, y=AllPC.2, color=AllPC.3)) + geom_text(label=rownames(All.scores))
sp+scale_color_gradient(low="red", high="darkgreen")
write.csv(All.scores, file = "NonphyloScores.csv")
#########################################################
############### Plotting #########
fviz_eig(pca)
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )    # Avoid text overlapping
biplot.phyl.pca()
biplot(princomp(log.WPTCS, cor = TRUE, scores = TRUE))

var <- get_pca_var(pca)
var
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)
# Coordinates of variables
head(var$coord, 4)

fviz_pca_var(pca, col.var = "blue")

fviz_pca_var(pca, labelsize = 5, repel = TRUE) 

  fviz_pca_var(pca, labelsize = 5, repel = TRUE)+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=96) +
  theme_minimal()
############################################################
  
#######################################################################
####################    Start here    #################################
#######################################################################
setwd("C:/Users/bfsco/Desktop/Masters Research/Color_analysis")

#### PGLS for-loop #### 
library(ape)
library(phytools)
library(Rcpp)
library(geiger)
library(RPANDA)
library(caper)
library(MCMCglmm)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(MuMIn)
Card.tree <- read.tree("Calibrated.Cardinal.tree")
row.names.remove <- c("Amaurospiza_concolor_aequatorialis",
                      "Chlorothraupis_carmioli_lutescens",
                      "Cyanoloxia_rothschildii",
                      "Granatellus_paraensis","Habia_fuscicauda_willisi",
                      "Habia_rubica_bahiae","Habia_rubica_peruviana",
                      "Pheucticus_aureoventris_uropygialis","Piranga_lutea_LSU_B5400",
                      "Pheucticus_chrysopeplus_aurantiacus")
Card.tree <- drop.tip(Card.tree,row.names.remove)  
row.names.remove2 <- c("Cyanoloxia_cyanoides")
Card.tree2 <- drop.tip(Card.tree,row.names.remove2)



#df <- read.csv("NEWDATA.noCy.csv", header = TRUE)
df <- read.csv("SEM.data.AllPCA.csv", header = TRUE)
DATA <- df
row.names(DATA) = DATA$Species
PTC <- scale(DATA$Percent.Tree.cover)
df <- cbind(df,PTC)
#get rid of node labels which prevent comparative dataframe object 
Card.tree$node.label<-NULL

comp.data<-comparative.data(Card.tree, df, names.col = "Species", vcv=TRUE,vcv.dim = 3, warn.dropped=TRUE)
full <- comp.data$data
Card.tree2$node.label<-NULL
comp.data2<-comparative.data(Card.tree2, df, names.col = "Species", vcv=TRUE,vcv.dim = 3, warn.dropped=TRUE)
#full2 <- comp.data2$data
##############################################3
  
a<-pgls(fAllPC1~mAllPC1, data=comp.data2,lambda = "ML") #
summary(a)
a$aicc

a<-pgls(mAllPC1~ClimComp.1+LandComp.3+NIR.mean, data=comp.data2,lambda = "ML") #
summary(a)
a$aicc
#############################################
a<-pgls(fAllPC1~VisRed.stdDev+Migration+Mass, data=comp.data,lambda = "ML") #
summary(a)
a$aicc

Strata
ClimComp.2
EVI.mean



with(df, plot(mAllPC2~ClimComp.2, xlab = "ClimComp.2", ylab = "mAllPC2",
              main = "mAllPC2 vs \ ClimComp.2"))
abline(a)

sp<-ggplot(full, aes(x=mAllPC1, y=mAllPC2, color=Percent.Tree.cover)) +  geom_text(label=rownames(full))
sp+scale_color_gradient( high="green",low="red")+  theme_light()


sp2<-ggplot(full, aes(x=fAllPC1, y=VisRed.stdDev)) +  geom_point(aes(colour = HabitatA), size = 3)
sp2+theme_bw() 
###############################
#Comparing all vs phylo

#df <- read.csv("SEM.Final.test.csv", header = TRUE)
phylo <- read.csv("SEM.Final.csv", header = TRUE)
phylo <- phylo[c(1:7)]
All <- cbind(df,phylo)
row.names(All) = All$Species

sp<-ggplot(All, aes(x=fAllPC1, y=fPC1, color=Percent.Tree.cover)) +  geom_text(label=rownames(All))
sp+scale_color_gradient( high="green",low="red")+  theme_light()

tTest <- t.test(All$mAllPC1,All$mPC1)
tTest$p.value
$statistic
