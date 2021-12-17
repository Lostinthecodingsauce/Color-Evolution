#Detecting outleirs 

#To identify the outliers in a PGLS model, first you need to extract the phylogenetic residuals from the model:
  
  resa<- residuals(a, phylo = TRUE) #extracts phylogenetic residuals from the pgls model
  
#Next these residuals need to be standardized by dividing by their square root of their variance:
    
    resa<- resa/sqrt(var(resa))[1] #standardises residuals by sqrt of their variance
  
  
#Finally the residuals can be matched to the species names and the outlier idenitified (in this case Colobus_polykomos):
    
    rownames(resa)<-rownames(a$residuals) #matches the residuals up with the species names
  
  rownames(resa)[(abs(resa)>3)]#gives the names of the outliers
  
# You can then remove these species from the data and tree and redo the analysis. Note that you may need to continue removing species until there are no more outliers.
  
  Cards1_nooutliers<-comp.data[-which(abs(resa)>3),]
  modela.pgls_nooutliers<-pgls(mPC1~ClimComp.2, data=Cards1_nooutliers, lambda =  "ML") ### Best primate_nooutliers, lambda='ML')
  summary(modela.pgls_nooutliers)
  plot(modela.pgls_nooutliers)
  abline(modela.pgls_nooutliers)
  
  
  ###########PC2 

  resb<- residuals(b, phylo = TRUE) #extracts phylogenetic residuals from the pgls model
  

  resb<- resb/sqrt(var(resb))[1] #standardises residuals by sqrt of their variance
  
  rownames(resb)<-rownames(b$residuals) #matches the residuals up with the species names
  
  rownames(resb)[(abs(resb)>3)]#gives the names of the outliers
  
  Cards_nooutliers<-This[-which(abs(resb)>3),]
  model.pgls_nooutliers<-pgls(mPC2~fPC2, data=Cards_nooutliers, lambda =  "ML") ### Best primate_nooutliers, lambda='ML')
  summary(model.pgls_nooutliers)
  plot(model.pgls_nooutliers)
  abline(model.pgls_nooutliers)
  
  ########### Other variables #########
  
  resb<- residuals(mPC1mods, phylo = TRUE) #extracts phylogenetic residuals from the pgls model
  
  
  resb<- resb/sqrt(var(resb))[1] #standardises residuals by sqrt of their variance
  
  rownames(resb)<-rownames(mPC1mods$residuals) #matches the residuals up with the species names
  
  rownames(resb)[(abs(resb)>3)]#gives the names of the outliers
  
  Cards_nooutliers<-df[-which(abs(resb)>3),]
  Cards_nooutliers<-comparative.data(Card.tree, Cards_nooutliers, names.col = "phylo", vcv=TRUE, warn.dropped=TRUE)
  
  
  model.pgls_nooutliers<-pgls(mPC1~Migration, data=Cards_nooutliers, lambda =  "ML") ### Best primate_nooutliers, lambda='ML')
  summary(model.pgls_nooutliers)
  plot(model.pgls_nooutliers)
  abline(model.pgls_nooutliers)