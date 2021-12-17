jetzDivRates <- function(tree) {
  
  spRate <- function(sp, tree) {
    #get branch lengths from root to tip
    edges <- vector()
    daughterFode <- match(sp, tree$tip.label)
    while (daughterFode != (length(tree$tip.label) + 1)) {
      parentFode <- tree$edge[which(tree$edge[,2] == daughterFode), 1]
      edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentFode & tree$edge[,2] == daughterFode)])
      daughterFode <- parentFode
    }
    
    res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
    res <- res ^ (-1)
    
    return(res)
  }
  
  rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
  names(rates) <- tree$tip.label
  
  return(rates)
  
}
