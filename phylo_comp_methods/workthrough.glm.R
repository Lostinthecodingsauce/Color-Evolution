
library(ape)
library(MCMCglmm)

########################################################################

setwd("C:/Users/bfsco/Desktop/Masters Research/Plumage/analysis/WPTSC/MT_tree_analysis")
#tip <- c("Cyanocompsa_cyanoides") #"Cyanocompsa_cyanoides"
FullCard.tree <- read.nexus("FullMCC.tree.nexus")
Card.tree <- read.nexus("FullMCC.tree.nexus")

All <- read.csv("new.csv", fileEncoding = 'UTF-8-BOM')

#####
# Set up ML PGLS Framework
##### 
library(phytools)
library(ape)
library(Rcpp)
library(geiger)
library(RPANDA)
library(caper)

data <- read.csv("new.csv", header = TRUE, fileEncoding = 'UTF-8-BOM') 
#Here I changed the names. Basically removed special characters, spaces, and renamed a couple columns for clarity
#I did this in excel because it was easier. 
names(data) <- c('X','fPC1','fPC2','fPC3','mPC1','mPC2','mPC3','Percent.Tree.cover','Habitat','Strata','StrataA','Forest.Dependency','Migration','MaleBrill','FemaleBrill',
                                     'NDVI.mean','NDVI.stdDev','NDVI.Var',
                                           'NIR.mean', 'NIR.stdDev','NIR.Var', 'GPP.mean','GPP.stdDev','GPP.Var') 



### Change to appropriate scale for spatial data. 
df.toscaled <- data[c(16:24)]*0.0001

temp <- subset(data[c(1:15)])
data <- cbind.data.frame(temp,df.toscaled) #df.toscale




#Analysis 1: GLMS for NDVI
#df <- df[,c(16:21)]

library(usdm)

df_cor <- data[c(16:24)]

#lets do some tests for colinnearity
cor(df_cor)
vif(df_cor)



#cor shows that Humidity and temperature are highly correlated (r = -.869)
#vif confirms this value (4 is the standard cutoff) and shows that Temperature has the most contribution to colinearity


#So now that we know our dataset a bit better, lets actually run some GLMS!
library(stats)
library(tidyverse)
library(MuMIn)

#get df without species level comparisons
df_m1 <-data[,c(5,8:13,16:24)]


df_m2 <-data[,c(1,6,9:13,16:24)]
df_f1 <-data[,c(1,2,9:13,16:24)]
df_f2 <-data[,c(1,3,9:13,16:24)]
#########################################

models <- list()
PCscores <- paste("data", 2:6, sep='')
ivnames <- paste("IV", 1:n, sep='') ## for some value of n

for (y in dvnames){
  form <- formula(paste(y, "~", ivnames))
  models[[y]] <- glm(form, data=dd, family='poisson') 
}
#To examine the results, just wrap your usual functions in a lapply, like this:
  
  lapply(models, summary) ## summarize each model










########################################

#this is a function that does iterative GLMs given a vector of formulas (i.e. Temp + Humidity, Altitude, Humidity*Temp+Altitude etc.)

multi_Glm <- function(data,respo = as.string(),vars){
  out <- list()
  for (i in 1:length(vars)){
    form = paste(respo, "~ ",vars[i])
  out[[i]] <- glm(form, data= data, family = "poisson") #%>%
  #anova(.,test = "Chisq")
  }
  names(out) <- vars
return(out)
}


#these lines will make every iteration of variables from all column names

#get column names of variables you want to test
vars <- colnames(df_m1[,5,8:24])


#Make every iteration using + as the seperator
formulas <- unlist(lapply(1:length(vars), function(i) combn(vars,i,paste,collapse = " + ")))

#do the GLMs using that custom function. You can change the response variable by changing the respo field to the #response column
a2<-multi_pgls(data=df_m1, respo = "mPC1", vars = formulas)


m1_glm <- multi_Glm(data = df_m1, respo = "Abundance", vars = formulas)
richness_glm  <- multi_Glm(data = df_sub, respo = "Richness", vars = formulas)
#here we make an AIC table to view all the formulas and find the best fit via AIC
#I also added BIC for better power and selecting models

aic_table_abundances <- data.frame("ID" = as.integer(),"form" = as.character(), "AICc" = as.integer(), "bic"= as.integer(), "deviance" = as.numeric())

for(i in 1:length(abundance_glm)){
  aic_table_abundances[i,1] <- i
aic_table_abundances[i,2] <- abundance_glm[[i]]$formula
aic_table_abundances[i,3] <- AICc(abundance_glm[[i]])
aic_table_abundances[i,4] <- BIC(abundance_glm[[i]])
aic_table_abundances[i,5] <- (1-(abundance_glm[[i]]$null.deviance/abundance_glm[[i]]$deviance))
}

aic_table_richness <- data.frame("ID" = as.integer(),"form" = as.character(), "AICc" = as.integer(), "bic"= as.integer(), "deviance" = as.numeric(),"delta_AICc" = as.numeric(), "delta_BIC" = as.numeric())

for(i in 1:length(richness_glm)){
  aic_table_richness[i,1] <- i
aic_table_richness[i,2] <- richness_glm[[i]]$formula
aic_table_richness[i,3] <- AICc(richness_glm[[i]])
aic_table_richness[i,4] <- BIC(richness_glm[[i]])
aic_table_richness[i,5] <- (1-(richness_glm[[i]]$null.deviance/richness_glm[[i]]$deviance))
}

#These will just give us the Deltas for each BIC and AICc
for(i in 1:nrow(aic_table_richness)){
  aic_table_richness[i,'delta_AICc'] <-aic_table_richness[i,3]- min(aic_table_richness$AICc)
   aic_table_richness[i,'delta_BIC'] <-aic_table_richness[i,4]- min(aic_table_richness$bic)
}

for(i in 1:nrow(aic_table_abundances)){
  aic_table_abundances[i,'delta_AICc'] <-aic_table_abundances[i,3]- min(aic_table_abundances$AICc)
   aic_table_abundances[i,'delta_BIC'] <-aic_table_abundances[i,4]- min(aic_table_abundances$bic)
}


```

```{r eval=FALSE, include=FALSE}
setwd('C:/Users/Admin/OneDrive - San Diego State University (SDSU.EDU)/Research/Personal Projects/EcuadorFrogs/EcuadorFrogs')

write.csv(aic_table_abundances, file = "Outputs/GLM_Abundance.csv")
write.csv(aic_table_richness, file = "Outputs/GLM_Richness.csv")
```

Here are the results from these GLMS. I included a metric of % deviance that is 1-null deviance/residual deviance. It gives us another idea of model fit. I sorted by BIC as that has an affinity for simpler models.


```{r include=FALSE}
abundance_top <- aic_table_abundances[order(aic_table_abundances$bic),]
richness_top <- aic_table_richness[order(aic_table_richness$bic),]

```

```{r}
richness_top[1:15,2:7]
abundance_top[1:15,2:7]
```
From the table we can see pretty clearly that Altitude keeps having a super large effect on the data for both abundances and richness. Same with SL. 

```{r eval=FALSE, include=FALSE}

summary(abundance_glm[[413]])

summary(richness_glm[[23]])


```





```{r eval=FALSE, include=FALSE}
#looking through we now do some summary state

lapply(test,anova,test = "Chisq")

summary(test[[413]])
anova(test[[413]],test = "Chisq")
lapply(test,summary)

#here im trying to get the values you got originally. 
t_glm <- glm(Abundance ~ Altitude_m + H + T + Rain_mm + H20_dist_m, data = df_sub, family = "poisson")
anova(t_glm, test = "Chisq")
summary(t_glm)

anova(test[[1]],test[[2]],test[[413]], test = "Chisq")


t_glm <- glm(Richness ~ CC + SL + HL + LL, data = df_sub, family = "poisson")
t_glm
summary(t_glm)
anova(t_glm,test = "Chisq")
plot(t_glm)
coefficients(t_glm)
summary(t_glm)




library(stargazer)

stargazer(aic_table_abundances[1:10,],digits = 2, 
          type = "html", 
          out = "Outputs/GLM_table.doc")
tab_model(aic_table_abundances)

model.sel(abundance_glm)
```


SAC

```{r}
library(vegan)
library(tidyverse)
library(plyr)

names(df)
#need to get the ataset with just abundances
df_sac <- df[,c(1,2,16:55)]
df_sac[is.na(df_sac)] <- 0

#next lets collapse based on Transect
df_sac <- aggregate(. ~Site+Transect, sum, data = df_sac)

sac <-specaccum(df_sac[,3:42], "random")

sac

plot(sac, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

df_sac_rando <- specaccum(df_sac[,3:42], "random")
plot(df_sac_rando,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Transects", ylab ="# of species found", main = "Species Accumulation Curve")
boxplot(df_sac_rando, col="yellow", add=TRUE, pch="+")

sites <- factor(df_sac$Site)
setwd('C:/Users/Admin/OneDrive - San Diego State University (SDSU.EDU)/Research/Personal Projects/EcuadorFrogs/EcuadorFrogs')
plot(df_sac_rando,  xlab = "Transects", ylab ="# of species found",  main = paste("Species Accumulation Curves per site"), xlim = c(0,12), ylim = c(0,30))

for(i in 1:3){
  df_sac_ <- specaccum(df_sac[which(df_sac$Site == sites[i]),3:42], "random")
plot(df_sac_, ci.lty=0, col= i, add = TRUE, lwd = 3)
write.csv(df_sac_,file = paste("Outputs/sac_site",i,".csv",sep = ""))
#boxplot(df_sac_, col="yellow", add=TRUE, pch="+")
}
text(locator(), labels = c("S1","S2","S3","Combined"))

plot(df_sac_rando)
plot(df_sac_, add = TRUE)
df_sac_ <- specaccum(df_sac[which(df_sac$Site = sites[1]),3:42], "random")

sites[i]
plot(df_sac_,ci.lty = 0, col = 3)
#Correlleogram

df_c <- df[,4:12]
names(df_c) <- c("CC","SL","HL","LL","T","H","R","A","D")


corrplot::corrplot(cor(df_c), type = "lower", order = "hclust", 
         tl.col = "black", tl.srt = 45, method = "color")
corrplot::corrplot.mixed(cor(df_c), lower.col = "black", number.cex = .7, upper = "color")
df_sac_ <- specaccum(df_sac[which(df_sac$Site == sites[3]),3:42], "random")
df_sac_

df[2,"Site"]



install.packages("visreg")
library(visreg)
library(ggplot2)
abundance_glm[100]

gg <- visreg::visreg(glm(Richness~Altitude_m+SL, data = df),
                     #"Altitude_m", 
                     overlay = FALSE, gg= TRUE, ylab = "Taxa Richness", xlab = "Altitude (m)") 


gg <- ggplot()+
  geom_point(data = df, aes(Altitude_m, Abundance, color = Site))+
  ylab("Abundance") + xlab("Altitude in meters")+
  geom_smooth(method = lm, aes(Altitude_m,Abundance), color = "dark grey", data = df)
gg


setwd('C:/Users/Admin/OneDrive - San Diego State University (SDSU.EDU)/Research/Personal Projects/EcuadorFrogs/EcuadorFrogs')
ggsave("Outputs/Alt_Abund_Color_Site.png", width =10, height = 7)


```
```{r}
library(fossil)
#Chao 1 gives an estimate of species richness given abundance data
#Lets do it for each site then

#divide the cleaned dataset by site. Turning NA into 0
#subseted datasets only include the abundance counts per species
site1 <- df[which(df$Site == "S1"),16:55]
site1[is.na(site1)] <- 0

site2 <- df[which(df$Site == "S2"),16:55]
site2[is.na(site2)] <- 0

site3 <- df[which(df$Site == "S3"),16:55]
site3[is.na(site3)] <- 0

#Run the Chao1 index across each site. Uses colsums to get the abundance per species
#uses the fossil library to do the calculations
#Values are commented beside each sites calculation
chao_s1 <- chao1(colSums(site1)) # 25.64
chao_s2 <- chao1(colSums(site2)) #26.08
chao_s3 <- chao1(colSums(site3)) #24.5

```

