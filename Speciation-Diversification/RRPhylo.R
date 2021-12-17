RRphylo <-
  function(tree,y,f=Ntip(tree)/10)
    ##### tree does not accept numbers as node labels ####
{
  require(ape)
  require(phytools)
  require(geiger)
  require(stats4)
  
  if (is.binary.tree(tree)) t <- tree else t <- multi2di(tree)
  
  internals <- unique(c(t$edge[, 1], t$edge[, 2][which(t$edge[, 
                                                              2] > Ntip(t))]))
  tippa <- list()
  for (i in 1:length(internals)) {
    tippas <- tips(t, internals[i])
    dato <- data.frame(rep(internals[i], length(tippas)), 
                       tippas)
    colnames(dato)[1] <- "node"
    tippa[[i]] <- dato
  }
  Tstr <- do.call(rbind, tippa)
  L <- matrix(nrow = Ntip(t), ncol = length(t$edge.length) + 
                1)
  rownames(L) <- t$tip.label
  colnames(L) <- c(internals, t$tip.label)
  edged <- data.frame(t$edge, t$edge.length)
  order <- edged[, 2][which(edged[, 2] < Ntip(t) + 1)]
  labs <- t$tip.label[order]
  edged[, 2][which(edged[, 2] < Ntip(t) + 1)] <- labs
  tip.path <- list()
  for (i in 1:dim(L)[1]) {
    tip.path[[i]] <- c(Tstr[, 1][which(Tstr[, 2] %in% rownames(L)[i])], 
                       which(t$tip.label == rownames(L)[i]))
  }
  for (j in 1:length(tip.path)) {
    a <- list()
    for (i in 2:length(tip.path[[j]]) - 1) {
      a[[i]] <- tip.path[[j]][c(i, i + 1)]
    }
    b <- do.call(rbind, a)
    b[which(b[, 2] < Ntip(t) + 1), 2] <- t$tip.label[b[which(b[, 
                                                               2] < Ntip(t) + 1), 2]]
    L.match <- b[, 2]
    br.len <- array()
    for (k in 1:dim(b)[1]) {
      br.len[k] <- edged[b[k, 1] == edged[, 1] & b[k, 2] == 
                           edged[, 2], ][, 3]
    }
    d <- data.frame(L.match, br.len)
    L[j, match(d[, 1], colnames(L))] <- d[, 2]
  }
  if (is.null(t$root)) L[, 1] <- 1 else L[, 1] <- t$root
  L[which(is.na(L))] <- 0
  edgedX <- data.frame(t$edge, t$edge.length)
  edged.1 <- edgedX[edgedX$X2 > Ntip(t), ]
  L1 <- matrix(ncol = length(internals), nrow = length(internals))
  colnames(L1) <- internals
  rownames(L1) <- internals
  node.path <- list()
  for (i in 1:length(internals)) {
    a <- getDescendants(t, internals[i])
    a <- a[a > Ntip(t)]
    node.path[[i]] <- data.frame(rep(internals[i], length(a)), 
                                 a)
  }
  node.path <- do.call(rbind, node.path)
  pathN <- list()
  for (i in 1:length(edged.1[, 2])) {
    pathN[[i]] <- c(node.path[which(node.path[, 2] == edged.1[i, 
                                                              2]), 1], edged.1[i, 2])
  }
  for (j in 1:length(pathN)) {
    a <- list()
    for (i in 2:length(pathN[[j]]) - 1) {
      a[[i]] <- pathN[[j]][c(i, i + 1)]
    }
    b <- do.call(rbind, a)
    L1.match <- b[, 2]
    br.len <- array()
    for (k in 1:dim(b)[1]) {
      br.len[k] <- edged.1[b[k, 1] == edged.1[, 1] & b[k, 
                                                       2] == edged.1[, 2], ][, 3]
    }
    d <- data.frame(L1.match, br.len)
    L1[j, match(d[, 1], colnames(L1))] <- d[, 2]
  }
  if (is.null(t$root)) L1[, 1] <- 1 else L1[, 1] <- t$root
  L1[which(is.na(L1))] <- 0
  internals <- unique(c(t$edge[, 1], t$edge[, 2][which(t$edge[, 
                                                              2] > Ntip(t))]))
  edged <- data.frame(t$edge, t$edge.length)
  optL <- function(lambda) {
    y <- scale(y)
    betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*% 
                t(L)) %*% as.matrix(y)
    aceRR <- L1 %*% betas[1:Nnode(t), ]
    y.hat <- L %*% betas
    Rvar <- array()
    for (i in 1:Ntip(t)) {
      ace.tip <- betas[match(names(which(L[i, ] != 0)), 
                             rownames(betas)), ]
      mat = as.matrix(dist(ace.tip))
      Rvar[i] <- sum(mat[row(mat) == col(mat) + 1])
    }
    abs(1 - (var(Rvar) + (mean(as.matrix(y))/mean(y.hat))))
  }
  lambda <- coef(mle(optL, start = list(lambda = 10), method = "L-BFGS-B"))
  betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*% 
              t(L)) %*% as.matrix(y)
  aceRR <- L1 %*% betas[1:Nnode(t), ]
  y.hat <- L %*% betas
  rates <- betas
  if (length(y) > Ntip(t)) {
    rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
    rates <- as.matrix(rates)
  } else { 
    rates <- rates
  }
  
  
  res <- list(t, L, L1, rates, aceRR, y.hat, lambda)
  names(res) <- c("tree", "tip.path", "node.path", "rates",  "aces", "y.estimates", "lambda")
  return(res)
}







search.shift <-function(RR,node=NULL,test.single=c("yes","no"),nrep=1000,status.type=c("clade","sparse"),state=NULL,auto.recognize=c("yes","no"),covariate=c("FALSE","TRUE"),cov=NULL,f=round(Ntip(tree)/10))
{
  require(phytools)
  require(geiger)
  require(scales)
  tree <- RR$tree
  rates <- RR$rates
  
  
  if (dim(rates)[2] > 1) {
    multi.rates <- apply(rates, 1, function(x) sqrt(sum(x^2)))
    rates <- as.matrix(multi.rates)
  } else {
    rates <- rates
  }
  
  match.arg(covariate)
  if(covariate=="FALSE"){
    rates<-rates
  }else{
    if(length(which(rates=="0"))>0){ 
      which(rates=="0")->zeroes
      log(abs(rates))->R
      R[-zeroes]->R
      
      if(length(cov)<length(R))
      {
        log(abs(c(RR$ace,cov)))->Y
        Y[-zeroes]->Y
        
      } else {
        cov->Y
        Y[-zeroes]->Y
      }
      residuals(lm(R~Y))->res
      which(rates!="0")->factOut
      
      rates[factOut]<-res
      rates[zeroes]<-0
    } else {
      log(abs(rates))->R
      log(abs(c(RR$ace,cov)))->Y
      residuals(lm(R~Y))->res
      as.matrix(res)->rates
    }
    
  }
  
  
  
  match.arg(status.type)
  if (status.type == "clade") {
    match.arg(auto.recognize)
    if (auto.recognize == "yes") {
      ST <- subtrees(tree)
      len <- array()
      for (i in 1:length(ST)) {
        len[i] <- Ntip(ST[[i]])
      }
      st <- ST[which(len < (Ntip(tree)/2) & len > round(f))]
      node <- sapply(st, function(x) getMRCA(tree, x$tip.label))
      names(st) <- node
      leaf2N.diff <- array()
      p.single <- array()
      for (j in 1:length(node)) {
        Cbranch <- getDescendants(tree, node[j])
        Ctips <- tips(tree, node[j])
        Cleaf <- c(Cbranch, Ctips)
        leaf.rates <- rates[match(Cleaf, rownames(rates)), 
                            ]
        leaf.rates <- na.omit(leaf.rates)
        NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
        leafR <- mean(abs(leaf.rates))
        NCR <- mean(abs(NCrates))
        leaf2N.diff[j] <- leafR - NCR
        NC <- length(rates) - length(leaf.rates)
        C <- length(leaf.rates)
        ran.diffM <- array()
        for (i in 1:nrep) {
          ran.diffM[i] <- mean(sample(abs(rates), C)) - mean(sample(abs(rates), 
                                                                    NC))
        }
        p.single[j] <- rank(c(leaf2N.diff[j], ran.diffM[1:(nrep - 
                                                             1)]))[1]/nrep
      }
      names(leaf2N.diff) <- node
      names(p.single) <- node
      if (length(p.single[p.single >= 0.975 | p.single <= 
                          0.025]) < 2) {
        p.single <- p.single[c(which.max(p.single), which.min(p.single))]
        leaf2N.diff <- leaf2N.diff[match(names(p.single), 
                                         names(leaf2N.diff))]
        p.init <- p.single
        l2N.init <- leaf2N.diff[match(names(p.init), 
                                      names(leaf2N.diff))]
      } else {
        p.single <- p.single[p.single >= 0.975 | p.single <= 0.025]
        
        leaf2N.diff <- leaf2N.diff[match(names(p.single), names(leaf2N.diff))]
        
        ups <- p.single[p.single > 0.975]
        dws <- p.single[p.single < 0.025]
        ups <- ups[na.omit(match(names(leaf2N.diff[order(leaf2N.diff, 
                                                         decreasing = FALSE)]), names(ups)))]
        dws <- dws[na.omit(match(names(leaf2N.diff[order(leaf2N.diff, 
                                                         decreasing = FALSE)]), names(dws)))]
        if (is.na(mean(dws))) {
          dws = Nnode(tree) * 2
        } else {
          
          
          s = 1 
          repeat
          { 
            d <- which(names(dws) %in% getDescendants(tree, names(dws)[s]))
            if (length(d) > 0) {
              leaf2N.diff[c(match(names(dws[d]),names(leaf2N.diff)),match(names(dws[s]),names(leaf2N.diff)))]->cla
              names(which.max(abs(leaf2N.diff[c(match(names(dws[d]),names(leaf2N.diff)),match(names(dws[s]),names(leaf2N.diff)))])))->IN
              dws[-match(names(cla[which(names(cla)!=IN)]),names(dws))]->dws
              s=1
            } else {
              dws <- dws
              s = s+1
              
            }
            
            if (s > length(dws))  break
            
          }
        }
        
        if (is.na(mean(ups))) {
          ups = Nnode(tree) * 2
        } else {
          
          z = 1
          repeat
          { 
            d <- which(names(ups) %in% getDescendants(tree, names(ups)[z]))
            if (length(d) > 0) {
              leaf2N.diff[c(match(names(ups[d]),names(leaf2N.diff)),match(names(ups[z]),names(leaf2N.diff)))]->cla
              names(which.max(abs(leaf2N.diff[c(match(names(ups[d]),names(leaf2N.diff)),match(names(ups[z]),names(leaf2N.diff)))])))->IN
              ups[-match(names(cla[which(names(cla)!=IN)]),names(ups))]->ups
              z=1
            } else {
              ups <- ups
              z = z+1
              
            }
            if (z > length(ups))  break
            
          }
        }
        
        
        p.init <- p.single
        l2N.init <- leaf2N.diff[match(names(p.init), 
                                      names(leaf2N.diff))]
        p.single <- p.single[which(names(p.single) %in% 
                                     names(c(ups, dws)))]
        
        leaf2N.diff <- leaf2N.diff[match(names(p.single), 
                                         names(leaf2N.diff))]
      }
      
      
      p.single[order(p.single)]->p.single
      pdf(file = "AR results for rate differences.pdf")
      if(Ntip(tree)>100) plot(tree, show.tip.label = FALSE) else plot(tree, cex=.8)
      xy <- list()
      for (w in 1:length(p.single)) {
        xy[[w]] <- unlist(sapply(get("last_plot.phylo", 
                                     envir = .PlotPhyloEnv), function(x) x[as.numeric(names(p.single)[w])]))[c(21, 
                                                                                                               22)]
      }
      
      
      c(rep("red",length(which(p.single<=0.025))),rep("royalblue",length(which(p.single>=0.975))))->p.col
      symbols(lapply(xy, "[[", 1), lapply(xy, "[[", 2), 
              circles = abs(leaf2N.diff[match(names(p.single),names(leaf2N.diff))])^0.5, inches = 0.25, 
              add = TRUE, bg = alpha(p.col, 0.5), fg = p.col)
      
      nodelabels(node = as.numeric(names(p.single)), adj = c(1.5, 
                                                             1), text = names(p.single), frame = "none", bg = "white", 
                 col = "purple")
      dev.off()
      p.shift = "AR mode, p.shift not implemented"
      return(list(`shift test` = p.shift, `all different clades` = p.init, 
                  `all clades rate difference` = l2N.init, `shift test single clades` = p.single, 
                  `average rate difference` = leaf2N.diff,`rates`=rates))
    } else {
      match.arg(test.single)
      if (test.single == "yes") {
        node = node
        Cbranch <- list()
        for (i in 1:length(node)) {
          Cbranch[[i]] <- getDescendants(tree, node[i])
        }
        Cbranch <- unlist(Cbranch)
        Cbranch <- Cbranch[-which(Cbranch < Ntip(tree))]
        Ctips <- list()
        for (i in 1:length(node)) {
          Ctips[[i]] <- tips(tree, node[i])
        }
        Ctips <- unlist(Ctips)
        Ctips <- unique(Ctips)
        Cbranch <- unique(Cbranch)
        Cleaf <- c(Cbranch, Ctips)
        leaf.rates <- rates[match(Cleaf, rownames(rates)), 
                            ]
        leaf.rates <- na.omit(leaf.rates)
        NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
        leafR <- mean(abs(leaf.rates))
        NCR <- mean(abs(NCrates))
        leaf2NC.diff <- leafR - NCR
        NC <- length(rates) - length(leaf.rates)
        C <- length(leaf.rates)
        ran.diffR <- array()
        for (i in 1:nrep) {
          ran.diffR[i] <- mean(sample(abs(rates), C)) - mean(sample(abs(rates), 
                                                                    NC))
        }
        p.shift <- rank(c(leaf2NC.diff, ran.diffR[1:(nrep - 
                                                       1)]))[1]/nrep
        par(mar = c(1, 1, 1, 1))
        par(mfrow = c(length(node) + 1, 1))
        hist(ran.diffR, xlab = "random differences", 
             main = "all clades", xlim = c(2.5 * range(ran.diffR)[1], 
                                           2.5 * range(ran.diffR)[2]))
        abline(v = leaf2NC.diff, col = "green", lwd = 3)
        leaf2N.diff <- array()
        p.single <- array()
        ran.diff <- list()
        for (i in 1:length(node)) {
          NOD <- node[-i]
          others <- list()
          mommies <- list()
          for (j in 1:length(NOD)) {
            others[[j]] <- tips(tree, NOD[j])
            mommies[[j]] <- getDescendants(tree, NOD[j])
          }
          others <- unlist(others)
          mommies <- unlist(mommies)
          mommies <- mommies[-which(mommies < Ntip(tree) + 
                                      1)]
          otmom <- c(mommies, others)
          Ctips <- tips(tree, node[i])
          Ctips <- unlist(Ctips)
          Cbranch <- getDescendants(tree, node[i])
          Cleaf <- c(Cbranch, Ctips)
          leaf.rates <- rates[match(Cleaf, rownames(rates)), 
                              ]
          leaf.rates <- na.omit(leaf.rates)
          NC <- rates[-c(which(rownames(rates) %in% names(leaf.rates)), 
                         which(rownames(rates) %in% otmom)), ]
          NR.r <- mean(abs(NC))
          leaf.r <- mean(abs(leaf.rates))
          leaf2N.diff[i] <- leaf.r - NR.r
          NC.l <- length(NC)
          leaf.l <- length(leaf.rates)
          tot.r <- abs(c(NC, leaf.rates))
          RAN.diff <- array()
          for (k in 1:nrep) {
            RAN.diff[k] <- mean(sample(tot.r, leaf.l)) - 
              mean(sample(tot.r, NC.l))
            ran.diff[[i]] <- RAN.diff
          }
          p.single[i] <- rank(c(leaf2N.diff[i], RAN.diff[1:(nrep - 
                                                              1)]))[1]/nrep
        }
        names(p.single) <- node
        names(leaf2N.diff) <- names(p.single)
        for (m in 1:length(node)) {
          hist(ran.diff[[m]], xlab = "random differences", 
               main = print(paste("Node", node[m], sep = " ")), 
               xlim = c(2.5 * range(ran.diff[[m]])[1], 2.5 * 
                          range(ran.diff[[m]])[2]))
          abline(v = leaf2N.diff[m], col = "blue", lwd = 3)
        }
        return(list(`shift test` = p.shift, 
                    `shift test single clades` = p.single, 
                    `average rate difference` = leaf2N.diff,`rates`=rates))
      } else {
        node = node
        Cbranch <- list()
        for (i in 1:length(node)) {
          Cbranch[[i]] <- getDescendants(tree, node[i])
        }
        Cbranch <- unlist(Cbranch)
        Cbranch <- Cbranch[-which(Cbranch < Ntip(tree))]
        Ctips <- list()
        for (i in 1:length(node)) {
          Ctips[[i]] <- tips(tree, node[i])
        }
        Ctips <- unlist(Ctips)
        Ctips <- unique(Ctips)
        Cbranch <- unique(Cbranch)
        Cleaf <- c(Cbranch, Ctips)
        leaf.rates <- rates[match(Cleaf, rownames(rates)), 
                            ]
        leaf.rates <- na.omit(leaf.rates)
        NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
        leafR <- mean(abs(leaf.rates))
        NCR <- mean(abs(NCrates))
        leaf2NC.diff <- leafR - NCR
        NC <- length(rates) - length(leaf.rates)
        C <- length(leaf.rates)
        ran.diffR <- array()
        for (i in 1:nrep) {
          ran.diffR[i] <- mean(sample(abs(rates), C)) - mean(sample(abs(rates), 
                                                                    NC))
        }
        p.shift <- rank(c(leaf2NC.diff, ran.diffR[1:(nrep - 
                                                       1)]))[1]/nrep
        hist(ran.diffR, xlab = "random differences", 
             main = "all clades", xlim = c(2.5 * range(ran.diffR)[1], 
                                           2.5 * range(ran.diffR)[2]))
        abline(v = leaf2NC.diff, col = "green", lwd = 3)
        return(list(`shift test` = p.shift, 
                    `average rate difference` = leaf2NC.diff,`rates`=rates))
      }
    }
  } else {
    frame <- data.frame(status = state, rate = rates[match(names(state), 
                                                           rownames(rates))])
    p.status.diff <- array()
    if (length(unique(state)) > 2) {
      status.diff <- apply(combn(tapply(abs(frame$rate), 
                                        frame$status, mean), 2), 2, diff)
      sta <- tapply(abs(frame$rate), frame$status, mean)
      sta <- sta[match(unique(state), names(sta))]
      w <- array()
      for (x in 1:length(sta)) w[x] <- sta[x] - mean(abs(frame[-which(frame$status == 
                                                                        names(sta)[x]), 2]))
      names(w) <- names(sta)
      status.diff <- c(status.diff, w)
      names(status.diff) <- c(apply(combn(levels(frame$status), 
                                          2), 2, function(x) paste(x[2], x[1], sep = "_")), 
                              names(sta))
      status.diffS <- matrix(ncol = length(status.diff), 
                             nrow = nrep)
      for (i in 1:nrep) {
        s.state <- frame$status
        s.ran <- sample(s.state)
        s.frame <- data.frame(s.ran, frame$rate)
        SD <- apply(combn(tapply(abs(s.frame$frame.rate), 
                                 s.frame$s.ran, mean), 2), 2, diff)
        sta <- tapply(abs(frame$rate), s.ran, mean)
        sta <- sta[match(unique(state), names(sta))]
        w <- array()
        for (x in 1:length(sta)) w[x] <- sta[x] - mean(abs(frame[-which(s.frame$s.ran == 
                                                                          names(sta)[x]), 2]))
        status.diffS[i, ] <- c(SD, w)
      }
      colnames(status.diffS) <- names(status.diff)
      par(mar = c(1, 1, 1, 1))
      par(mfrow = c(length(unique(state)), 1))
      idx <- match(unique(state), colnames(status.diffS))
      for (i in 1:length(idx)) {
        hist(status.diffS[, idx[i]], xlab = "random differences", 
             main = print(paste("rate difference per status", 
                                colnames(status.diffS)[idx[i]], sep = " ")), 
             xlim = c(min(status.diffS[, idx[i]]) - sd(status.diffS[, 
                                                                    idx[i]]), max(status.diffS[, idx[i]]) + sd(status.diffS[, 
                                                                                                                            idx[i]])))
        abline(v = status.diff[idx[i]], lwd = 3, col = "green")
      }
      for (i in 1:length(status.diff)) p.status.diff[i] <- rank(c(status.diff[i], 
                                                                  status.diffS[1:(nrep - 1), i]))[1]/nrep
      names(p.status.diff) <- names(status.diff)
      return(list(`rate difference` = status.diff, `p rate diff` = p.status.diff,`rates`=rates))
    } else {
      status.diff <- diff(tapply(abs(frame$rate), state, 
                                 mean))
      status.diffS <- array()
      for (i in 1:nrep) {
        s.state <- frame$status
        s.ran <- sample(s.state)
        s.frame <- data.frame(s.ran, frame$rate)
        s.frame[, 1] <- as.factor(s.frame[, 1])
        status.diffS[i] <- diff(tapply(abs(s.frame$frame.rate), 
                                       s.frame$s.ran, mean))
      }
      hist(status.diffS, xlab = "random differences", main = "rate difference per status", 
           xlim = c(min(status.diffS) * 2.5, max(status.diffS) * 
                      2.5))
      abline(v = status.diff, lwd = 3, col = "green")
      p.status.diff <- rank(c(status.diff, status.diffS[1:(nrep - 
                                                             1)]))[1]/nrep
      return(list(`rate difference` = status.diff, `p rate diff` = p.status.diff,`rates`=rates))
    }
  }
}
setBM <-
  function(tree,nY=1,s2=1,a=0,type=c("","brown","trend","drift"),tr=5,t.shift=.5,trend.type=c("linear","stepwise"))
  {
    require(lmtest)
    require(phytools)
    type <- match.arg(type)
    if (type == "") 
      stop("argument ‘type’ must be defined")
    switch(type, brown = {
      i = 1
      times <- diag(vcv(tree))
      yy <- list()
      while (length(yy) < nY) {
        y <- fastBM(tree, sig2 = s2, a = a)
        res <- bptest(y ~ times)[4]
        if (res > 0.05 & coefficients(summary(lm(y ~ times)))[8] > 
            0.05) {
          yy[[i]] <- y
        }
        yy <- Filter(Negate(is.null), yy)
        i = i + 1
      }
      yy <- do.call(rbind, yy)
      yy <- t(yy)
      if (nY == 1) {
        nam <- rownames(yy)
        yy <- array(yy)
        names(yy) <- nam
      }
    }, trend = {
      i = 1
      times <- diag(vcv(tree))
      yy <- list()
      while (length(yy) < nY) {
        y <- fastBM(tree, sig2 = s2, a = a)
        res <- bptest(y ~ times)[4]
        if (res > 0.05 & coefficients(summary(lm(y ~ times)))[8] > 
            0.05) {
          match.arg(trend.type)
          if (trend.type == "linear") {
            y <- diag(vcv(tree))/min(diag(vcv(tree))) * 
              y
          } else {
            y[match(names(which(times > t.shift * max(nodeHeights(tree)))), 
                    names(y))] <- y[match(names(which(times > 
                                                        t.shift * max(nodeHeights(tree)))), names(y))] * 
              tr
          }
          yy[[i]] <- y
        }
        yy <- Filter(Negate(is.null), yy)
        i = i + 1
      }
      yy <- do.call(rbind, yy)
      yy <- t(yy)
      if (nY == 1) {
        nam <- rownames(yy)
        yy <- array(yy)
        names(yy) <- nam
      }
    }, drift = {
      i = 1
      times <- diag(vcv(tree))
      yy <- list()
      while (length(yy) < nY) {
        y <- fastBM(tree, sig2 = s2, a = a)
        res <- bptest(y ~ times)[4]
        if (res > 0.05 & coefficients(summary(lm(y ~ times)))[8] < 
            0.001) {
          yy[[i]] <- y
        }
        yy <- Filter(Negate(is.null), yy)
        i = i + 1
      }
      yy <- do.call(rbind, yy)
      yy <- t(yy)
      if (nY == 1) {
        nam <- rownames(yy)
        yy <- array(yy)
        names(yy) <- nam
      }
    })
    return(yy)
  }
sizedsubtree <-function(tree,size=Ntip(tree)/10)
{
  require(ape)
  require(geiger)
  nod <- array()
  repeat {
    i = 1
    node <- sample(seq(Ntip(tree)+2, dim(tree$edge)[1]-1), 1)
    a <- length(tips(tree, node))
    i = i + 1
    if (a >= size & a <= (Ntip(tree)-1)/2) 
      nod[i] <- node
    else nod[i] <- 1
    if (nod[i] > 1) {
      break
    }
    if (is.na(sum(nod))) nod <- na.omit(nod) else nod <- nod
    
  }
  return(nod[2])
}

