#### Raw Color Values PGLS ####
### THis script is the same as the Phylo.linear.mixed models script except that I am using log-transformed raw traits rther than PCA

#########################################################
setwd("C:/Users/bfsco/Desktop/Masters Research/Plumage/analysis/WPTSC/MT_tree_analysis")

##### ML PGLS
library(phytools)
library(ape)
library(Rcpp)
library(geiger)
library(RPANDA)
library(caper)

data <- read.csv("RawTraits_NDVI.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') # Non-scaled spatial data
FullCard.tree <- read.nexus("FullMCC.tree.nexus")


# Transform Spatial data to have correct scale 
df.scale <- data[c(31:36)]*0.0001
temp <- subset(data[c(1:29)])
data <- cbind.data.frame(temp,df.scale)

comp.data<-comparative.data(FullCard.tree, data, names.col = "Species", vcv=TRUE, warn.dropped=TRUE)

###########
#NDVI.mean	NDVI.stdDev	NDVI.min	NDVI.max	NDVI.skew	
#EVI.mean	EVI.stdDev	EVI.min	EVI.max	EVI.skew

#Percent.Tree.cover
mVol<-pgls(M_Volume~Percent.Tree.cover, data=comp.data,lambda = "ML") ### 
fVol<-pgls(F_Volume~Percent.Tree.cover, data=comp.data,lambda = "ML") ## 

mAvgHueDis<-pgls(M_AvgHueDisp~Percent.Tree.cover, data=comp.data,lambda = "ML")## 
fAvgHueDis<-pgls(F_AvgHueDisp~Percent.Tree.cover, data=comp.data,lambda = "ML")## 
mMaxHueDis<-pgls(M_MaxHueDisp~Percent.Tree.cover, data=comp.data,lambda = "ML")## 
fMaxHueDis<-pgls(F_MaxHueDisp~Percent.Tree.cover, data=comp.data,lambda = "ML")## 
mVarHueDis<-pgls(M_VarHueDisp~Percent.Tree.cover, data=comp.data,lambda = "ML") ### 
fVarHueDis<-pgls(F_VarHueDisp~Percent.Tree.cover, data=comp.data,lambda = "ML") ### 

mAvgSpan<-pgls(M_AvgSpan~Percent.Tree.cover, data=comp.data,lambda = "ML") ### 
fAvgSpan<-pgls(F_AvgSpan~Percent.Tree.cover, data=comp.data,lambda = "ML") ## 
mMaxSpan<-pgls(M_MaxSpan~Percent.Tree.cover, data=comp.data,lambda = "ML") ### 
fMaxSpan<-pgls(F_MaxSpan~Percent.Tree.cover, data=comp.data,lambda = "ML") ## 
mVarSpan<-pgls(M_VarSpan~Percent.Tree.cover, data=comp.data,lambda = "ML") ### 
fVarSpan<-pgls(F_VarSpan~Percent.Tree.cover, data=comp.data,lambda = "ML") ## 

mChroma<-pgls(M_AvgChroma~Percent.Tree.cover, data=comp.data,lambda = "ML")## 
fChroma<-pgls(F_AvgChroma~Percent.Tree.cover, data=comp.data,lambda = "ML")## 

mBrill<-pgls(M_AvgBrill~Percent.Tree.cover, data=comp.data,lambda = "ML")## 
fBrill<-pgls(F_AvgBrill~Percent.Tree.cover, data=comp.data,lambda = "ML")## 


summary(mVol) 
summary(fVol)

summary(mAvgHueDis)
summary(fAvgHueDis)

summary(mMaxHueDis)
summary(fMaxHueDis)

summary(mVarHueDis)
summary(fVarHueDis)

summary(mAvgSpan)
summary(fAvgSpan)

summary(mMaxSpan)
summary(fMaxSpan) #yes

summary(mVarSpan)
summary(fVarSpan) #yes 


summary(mChroma)
summary(fChroma)


summary(mBrill) #sig 
summary(fBrill) #sig


#NDVI.mean
mVol<-pgls(M_Volume~NDVI.mean, data=comp.data,lambda = "ML") ### 
fVol<-pgls(F_Volume~NDVI.mean, data=comp.data,lambda = "ML") ## 

mAvgHueDis<-pgls(M_AvgHueDisp~NDVI.mean, data=comp.data,lambda = "ML")## 
fAvgHueDis<-pgls(F_AvgHueDisp~NDVI.mean, data=comp.data,lambda = "ML")## 
mMaxHueDis<-pgls(M_MaxHueDisp~NDVI.mean, data=comp.data,lambda = "ML")## 
fMaxHueDis<-pgls(F_MaxHueDisp~NDVI.mean, data=comp.data,lambda = "ML")## 
mVarHueDis<-pgls(M_VarHueDisp~NDVI.mean, data=comp.data,lambda = "ML") ### 
fVarHueDis<-pgls(F_VarHueDisp~NDVI.mean, data=comp.data,lambda = "ML") ### 

mAvgSpan<-pgls(M_AvgSpan~NDVI.mean, data=comp.data,lambda = "ML") ### 
fAvgSpan<-pgls(F_AvgSpan~NDVI.mean, data=comp.data,lambda = "ML") ## 
mMaxSpan<-pgls(M_MaxSpan~NDVI.mean, data=comp.data,lambda = "ML") ### 
fMaxSpan<-pgls(F_MaxSpan~NDVI.mean, data=comp.data,lambda = "ML") ## 
mVarSpan<-pgls(M_VarSpan~NDVI.mean, data=comp.data,lambda = "ML") ### 
fVarSpan<-pgls(F_VarSpan~NDVI.mean, data=comp.data,lambda = "ML") ## 

mChroma<-pgls(M_AvgChroma~NDVI.mean, data=comp.data,lambda = "ML")## 
fChroma<-pgls(F_AvgChroma~NDVI.mean, data=comp.data,lambda = "ML")## 

mBrill<-pgls(M_AvgBrill~NDVI.mean, data=comp.data,lambda = "ML")## 
fBrill<-pgls(F_AvgBrill~NDVI.mean, data=comp.data,lambda = "ML")## 


summary(mVol) # yes
summary(fVol)

summary(mAvgHueDis)
summary(fAvgHueDis)

summary(mMaxHueDis)
summary(fMaxHueDis)

summary(mVarHueDis)
summary(fVarHueDis)

summary(mAvgSpan)
summary(fAvgSpan)

summary(mMaxSpan)
summary(fMaxSpan) #yes

summary(mVarSpan)
summary(mVarSpan)


summary(mChroma)
summary(fChroma)


summary(mBrill) #sig 
summary(fBrill) #sig

#######################
# VisRed old
####################
#VisRed.mean Old
mVol<-pgls(M_Volume~VisRed.mOld, data=comp.data,lambda = "ML") ### 
fVol<-pgls(F_Volume~VisRed.mOld, data=comp.data,lambda = "ML") ## 

mAvgHueDis<-pgls(M_AvgHueDisp~VisRed.mOld, data=comp.data,lambda = "ML")## 
fAvgHueDis<-pgls(F_AvgHueDisp~VisRed.mOld, data=comp.data,lambda = "ML")## 
mMaxHueDis<-pgls(M_MaxHueDisp~VisRed.mOld, data=comp.data,lambda = "ML")## 
fMaxHueDis<-pgls(F_MaxHueDisp~VisRed.mOld, data=comp.data,lambda = "ML")## 
mVarHueDis<-pgls(M_VarHueDisp~VisRed.mOld, data=comp.data,lambda = "ML") ### 
fVarHueDis<-pgls(F_VarHueDisp~VisRed.mOld, data=comp.data,lambda = "ML") ### 

mAvgSpan<-pgls(M_AvgSpan~VisRed.mOld, data=comp.data,lambda = "ML") ### 
fAvgSpan<-pgls(F_AvgSpan~VisRed.mOld, data=comp.data,lambda = "ML") ## 
mMaxSpan<-pgls(M_MaxSpan~VisRed.mOld, data=comp.data,lambda = "ML") ### 
fMaxSpan<-pgls(F_MaxSpan~VisRed.mOld, data=comp.data,lambda = "ML") ## 
mVarSpan<-pgls(M_VarSpan~VisRed.mOld, data=comp.data,lambda = "ML") ### 
fVarSpan<-pgls(F_VarSpan~VisRed.mOld, data=comp.data,lambda = "ML") ## 

mChroma<-pgls(M_AvgChroma~VisRed.mOld, data=comp.data,lambda = "ML")## 
fChroma<-pgls(F_AvgChroma~VisRed.mOld, data=comp.data,lambda = "ML")## 

mBrill<-pgls(M_AvgBrill~VisRed.mOld, data=comp.data,lambda = "ML")## 
fBrill<-pgls(F_AvgBrill~VisRed.mOld, data=comp.data,lambda = "ML")## 


summary(mVol) 
summary(fVol) # yes

summary(mAvgHueDis)
summary(fAvgHueDis)

summary(mMaxHueDis)
summary(fMaxHueDis)

summary(mVarHueDis)
summary(fVarHueDis) # yes

summary(mAvgSpan)
summary(fAvgSpan) # yes

summary(mMaxSpan)
summary(fMaxSpan) # yes

summary(mVarSpan)
summary(fVarSpan) # yes


summary(mChroma)
summary(fChroma)


summary(mBrill) #yes
summary(fBrill) #yes 


##################################3
 # Vis Red New
######################################
#VisRed.mean 
mVol<-pgls(M_Volume~VisRed.mean, data=comp.data,lambda = "ML") ### 
fVol<-pgls(F_Volume~VisRed.mean, data=comp.data,lambda = "ML") ## 

mAvgHueDis<-pgls(M_AvgHueDisp~VisRed.mean, data=comp.data,lambda = "ML")## 
fAvgHueDis<-pgls(F_AvgHueDisp~VisRed.mean, data=comp.data,lambda = "ML")## 
mMaxHueDis<-pgls(M_MaxHueDisp~VisRed.mean, data=comp.data,lambda = "ML")## 
fMaxHueDis<-pgls(F_MaxHueDisp~VisRed.mean, data=comp.data,lambda = "ML")## 
mVarHueDis<-pgls(M_VarHueDisp~VisRed.mean, data=comp.data,lambda = "ML") ### 
fVarHueDis<-pgls(F_VarHueDisp~VisRed.mean, data=comp.data,lambda = "ML") ### 

mAvgSpan<-pgls(M_AvgSpan~VisRed.mean, data=comp.data,lambda = "ML") ### 
fAvgSpan<-pgls(F_AvgSpan~VisRed.mean, data=comp.data,lambda = "ML") ## 
mMaxSpan<-pgls(M_MaxSpan~VisRed.mean, data=comp.data,lambda = "ML") ### 
fMaxSpan<-pgls(F_MaxSpan~VisRed.mean, data=comp.data,lambda = "ML") ## 
mVarSpan<-pgls(M_VarSpan~VisRed.mean, data=comp.data,lambda = "ML") ### 
fVarSpan<-pgls(F_VarSpan~VisRed.mean, data=comp.data,lambda = "ML") ## 

mChroma<-pgls(M_AvgChroma~VisRed.mean, data=comp.data,lambda = "ML")## 
fChroma<-pgls(F_AvgChroma~VisRed.mean, data=comp.data,lambda = "ML")## 

mBrill<-pgls(M_AvgBrill~VisRed.mean, data=comp.data,lambda = "ML")## 
fBrill<-pgls(F_AvgBrill~VisRed.mean, data=comp.data,lambda = "ML")## 


summary(mVol) 
summary(fVol)

summary(mAvgHueDis)
summary(fAvgHueDis)

summary(mMaxHueDis)
summary(fMaxHueDis)

summary(mVarHueDis)
summary(fVarHueDis)

summary(mAvgSpan)
summary(fAvgSpan) # yes

summary(mMaxSpan)
summary(fMaxSpan) ## yes

summary(mVarSpan)
summary(fVarSpan) ## yes


summary(mChroma)
summary(fChroma)


summary(mBrill) 
summary(fBrill) # yes 

####################### Plotting ##################
library(ggimage)
library(ggtree)
library(ggplot2)


sp2<-ggplot(data, aes(x=fMaxSpan, y= Percent.Tree.cover, color=VisRed.mean)) + geom_point(size = 3)
sp2#+scale_color_gradient(low="blue", high="red")





################# Standard DEvation ###########################
#NDVI.stdDev
mVol<-pgls(M_Volume~NDVI.stdDev, data=comp.data,lambda = "ML") ### 
fVol<-pgls(F_Volume~NDVI.stdDev, data=comp.data,lambda = "ML") ## 

mAvgHueDis<-pgls(M_AvgHueDisp~NDVI.stdDev, data=comp.data,lambda = "ML")## 
fAvgHueDis<-pgls(F_AvgHueDisp~NDVI.stdDev, data=comp.data,lambda = "ML")## 
mMaxHueDis<-pgls(M_MaxHueDisp~NDVI.stdDev, data=comp.data,lambda = "ML")## 
fMaxHueDis<-pgls(F_MaxHueDisp~NDVI.stdDev, data=comp.data,lambda = "ML")## 
mVarHueDis<-pgls(M_VarHueDisp~NDVI.stdDev, data=comp.data,lambda = "ML") ### 
fVarHueDis<-pgls(F_VarHueDisp~NDVI.stdDev, data=comp.data,lambda = "ML") ### 

mAvgSpan<-pgls(M_AvgSpan~NDVI.stdDev, data=comp.data,lambda = "ML") ### 
fAvgSpan<-pgls(F_AvgSpan~NDVI.stdDev, data=comp.data,lambda = "ML") ## 
mMaxSpan<-pgls(M_MaxSpan~NDVI.stdDev, data=comp.data,lambda = "ML") ### 
fMaxSpan<-pgls(F_MaxSpan~NDVI.stdDev, data=comp.data,lambda = "ML") ## 
mVarSpan<-pgls(M_VarSpan~NDVI.stdDev, data=comp.data,lambda = "ML") ### 
fVarSpan<-pgls(F_VarSpan~NDVI.stdDev, data=comp.data,lambda = "ML") ## 

mChroma<-pgls(M_AvgChroma~NDVI.stdDev, data=comp.data,lambda = "ML")## 
fChroma<-pgls(F_AvgChroma~NDVI.stdDev, data=comp.data,lambda = "ML")## 

mBrill<-pgls(M_AvgBrill~NDVI.stdDev, data=comp.data,lambda = "ML")## 
fBrill<-pgls(F_AvgBrill~NDVI.stdDev, data=comp.data,lambda = "ML")## 


summary(mVol) 
summary(fVol)

summary(mAvgHueDis)
summary(fAvgHueDis)

summary(mMaxHueDis)
summary(fMaxHueDis)

summary(mVarHueDis)
summary(fVarHueDis)

summary(mAvgSpan)
summary(fAvgSpan)

summary(mMaxSpan)
summary(fMaxSpan) 

summary(mVarSpan)
summary(fVarSpan)


summary(mChroma)
summary(fChroma)


summary(mBrill) 
summary(fBrill) # yes

#######################
# VisRed old
####################
#VisRed.stdDev Old
mVol<-pgls(M_Volume~VisRed.stdOld, data=comp.data,lambda = "ML") ### 
fVol<-pgls(F_Volume~VisRed.stdOld, data=comp.data,lambda = "ML") ## 

mAvgHueDis<-pgls(M_AvgHueDisp~VisRed.stdOld, data=comp.data,lambda = "ML")## 
fAvgHueDis<-pgls(F_AvgHueDisp~VisRed.stdOld, data=comp.data,lambda = "ML")## 
mMaxHueDis<-pgls(M_MaxHueDisp~VisRed.stdOld, data=comp.data,lambda = "ML")## 
fMaxHueDis<-pgls(F_MaxHueDisp~VisRed.stdOld, data=comp.data,lambda = "ML")## 
mVarHueDis<-pgls(M_VarHueDisp~VisRed.stdOld, data=comp.data,lambda = "ML") ### 
fVarHueDis<-pgls(F_VarHueDisp~VisRed.stdOld, data=comp.data,lambda = "ML") ### 

mAvgSpan<-pgls(M_AvgSpan~VisRed.stdOld, data=comp.data,lambda = "ML") ### 
fAvgSpan<-pgls(F_AvgSpan~VisRed.stdOld, data=comp.data,lambda = "ML") ## 
mMaxSpan<-pgls(M_MaxSpan~VisRed.stdOld, data=comp.data,lambda = "ML") ### 
fMaxSpan<-pgls(F_MaxSpan~VisRed.stdOld, data=comp.data,lambda = "ML") ## 
mVarSpan<-pgls(M_VarSpan~VisRed.stdOld, data=comp.data,lambda = "ML") ### 
fVarSpan<-pgls(F_VarSpan~VisRed.stdOld, data=comp.data,lambda = "ML") ## 

mChroma<-pgls(M_AvgChroma~VisRed.stdOld, data=comp.data,lambda = "ML")## 
fChroma<-pgls(F_AvgChroma~VisRed.stdOld, data=comp.data,lambda = "ML")## 

mBrill<-pgls(M_AvgBrill~VisRed.stdOld, data=comp.data,lambda = "ML")## 
fBrill<-pgls(F_AvgBrill~VisRed.stdOld, data=comp.data,lambda = "ML")## 


summary(mVol) 
summary(fVol) # yes

summary(mAvgHueDis)
summary(fAvgHueDis)

summary(mMaxHueDis)
summary(fMaxHueDis)

summary(mVarHueDis)
summary(fVarHueDis) # yes

summary(mAvgSpan)
summary(fAvgSpan) # yes

summary(mMaxSpan)
summary(fMaxSpan) # yes

summary(mVarSpan)
summary(fVarSpan) # yes


summary(mChroma)
summary(fChroma)


summary(mBrill) #yes
summary(fBrill) #yes 


##################################3
# Vis Red New
######################################
#VisRed.stdDev 
mVol<-pgls(M_Volume~VisRed.stdDev, data=comp.data,lambda = "ML") ### 
fVol<-pgls(F_Volume~VisRed.stdDev, data=comp.data,lambda = "ML") ## 

mAvgHueDis<-pgls(M_AvgHueDisp~VisRed.stdDev, data=comp.data,lambda = "ML")## 
fAvgHueDis<-pgls(F_AvgHueDisp~VisRed.stdDev, data=comp.data,lambda = "ML")## 
mMaxHueDis<-pgls(M_MaxHueDisp~VisRed.stdDev, data=comp.data,lambda = "ML")## 
fMaxHueDis<-pgls(F_MaxHueDisp~VisRed.stdDev, data=comp.data,lambda = "ML")## 
mVarHueDis<-pgls(M_VarHueDisp~VisRed.stdDev, data=comp.data,lambda = "ML") ### 
fVarHueDis<-pgls(F_VarHueDisp~VisRed.stdDev, data=comp.data,lambda = "ML") ### 

mAvgSpan<-pgls(M_AvgSpan~VisRed.stdDev, data=comp.data,lambda = "ML") ### 
fAvgSpan<-pgls(F_AvgSpan~VisRed.stdDev, data=comp.data,lambda = "ML") ## 
mMaxSpan<-pgls(M_MaxSpan~VisRed.stdDev, data=comp.data,lambda = "ML") ### 
fMaxSpan<-pgls(F_MaxSpan~VisRed.stdDev, data=comp.data,lambda = "ML") ## 
mVarSpan<-pgls(M_VarSpan~VisRed.stdDev, data=comp.data,lambda = "ML") ### 
fVarSpan<-pgls(F_VarSpan~VisRed.stdDev, data=comp.data,lambda = "ML") ## 

mChroma<-pgls(M_AvgChroma~VisRed.stdDev, data=comp.data,lambda = "ML")## 
fChroma<-pgls(F_AvgChroma~VisRed.stdDev, data=comp.data,lambda = "ML")## 

mBrill<-pgls(M_AvgBrill~VisRed.stdDev, data=comp.data,lambda = "ML")## 
fBrill<-pgls(F_AvgBrill~VisRed.stdDev, data=comp.data,lambda = "ML")## 


summary(mVol) 
summary(fVol)

summary(mAvgHueDis)
summary(fAvgHueDis)

summary(mMaxHueDis)
summary(fMaxHueDis)

summary(mVarHueDis)
summary(fVarHueDis)

summary(mAvgSpan)
summary(fAvgSpan) # yes

summary(mMaxSpan)
summary(fMaxSpan) ## yes

summary(mVarSpan)
summary(fVarSpan) ## yes


summary(mChroma)
summary(fChroma)


summary(mBrill) 
summary(fBrill) # yes 
