
plotRates<-function(RR, node, export.tiff=c(TRUE, FALSE),foldername){
  #require(phytools)
  
  RR$rates->FULLrates
  RR$tree->rr3
  getDescendants(rr3,node)->shift331
  
  c(FULLrates[match(shift331[shift331>Ntip(rr3)],rownames(FULLrates)),],FULLrates[match(rr3$tip.label[shift331[shift331<=Ntip(rr3)]],rownames(FULLrates)),])->shift.rates
  shift.rates[order(shift.rates,decreasing=TRUE)]->shift.rates
  
  
  
  if(export.tiff==TRUE)
  {
    tiff(paste(foldername,"rate_bars.tiff",sep="/"),
         width=2000,
         height=1500,
         pointsize=4,
         res=600,
         compression="lzw")
    par(mfrow=c(1,2))
    hist(log(abs(shift.rates)))->H1
    hist(log(abs(FULLrates[-match(names(shift.rates),rownames(FULLrates)),])))->H2
    log(abs(FULLrates))[(log(abs(FULLrates))!="-Inf")]->Xa
    par(mar=c(4,4,4,4))
    plot(H2,col=rgb(0,1,0,.75),xlim=c(range(Xa)[1]-diff(range(Xa))*.1,range(Xa)[2]+diff(range(Xa))*.1),ylim=c(0,max(H2$count)*1.5),xlab="log absolute rates",main="")
    plot(H1,col=rgb(0,0,1,1),add=TRUE)
    legend("topleft", c("back rates", "shift node rates"), col=c(rgb(0,1,0,.75),rgb(0,0,1,1)), box.lwd = 0,box.col = "white",bg = "white",lwd=10)
    par(las=2)
    par(mar=c(4,10,4,4))
    barplot(shift.rates,xlab="rates",horiz=TRUE,col="blue",names.arg=names(shift.rates), main="",border="red")
    abline(v=mean(FULLrates),col="red",lwd=3)
    
    dev.off()
    return(shift.rates)
  }else{
    
    par(mfrow=c(1,2))
    hist(log(abs(shift.rates)))->H1
    hist(log(abs(FULLrates[-match(names(shift.rates),rownames(FULLrates)),])))->H2
    log(abs(FULLrates))[(log(abs(FULLrates))!="-Inf")]->Xa
    par(mar=c(4,4,4,4))
    plot(H2,col=rgb(0,1,0,.75),xlim=c(range(Xa)[1]-diff(range(Xa))*.1,range(Xa)[2]+diff(range(Xa))*.1),ylim=c(0,max(H2$count)*1.3),xlab="log absolute rates",main="")
    plot(H1,col=rgb(0,0,1,1),add=TRUE)
    legend("topleft", c("back rates", "shift node rates"), col=c(rgb(0,1,0,.75),rgb(0,0,1,1)), box.lwd = 0,box.col = "white",bg = "white",lwd=10)
    par(las=2)
    par(mar=c(4,10,4,4))
    barplot(shift.rates,xlab="rates",horiz=TRUE,col="blue",names.arg=names(shift.rates), main="",border="red")
    abline(v=mean(FULLrates),col="red",lwd=3)
    
    return(shift.rates)
  }
}