# load tree
tree <- read.tree("Thraupidae_pruned_21_01_2019.tre")

# load phylo PC scores
phy.pc.scores <- readRDS("Thraupidae_pruned_bill_9_pPCs.rds")

#########################################
##Quantify morphological specialization##
#########################################

# position in morphospace is the square root of the sum of squares for each score (Ricklefs, 2012)
position.in.morphospace <- function(scores){
  scores <- scores^2
  pos <- scores %>% sum %>% sqrt
  return(pos)
}

# calculate distance from center of morphospace
phy.D.scores <- apply(phy.pc.scores, 1, position.in.morphospace)