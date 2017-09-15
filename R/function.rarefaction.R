################################################################################
#Toolbox 16S
#
#Convenience functions to perform rarefaction analyses
#
#2017-09-15
#sebastian.schmidt@embl.de
################################################################################

################################################################################
################################################################################
#Function: rarefaction_analysis
#
#Input: count table, (relative) rarefaction steps, iterations per step (optional), Hill q to perform
################################################################################
################################################################################
source("https://raw.githubusercontent.com/defleury/Toolbox_16S/master/R/function.alpha_diversity.R")
rarefaction <- function(count.table, steps=c(seq(0.01, 0.09, by=0.01), seq(0.1, 0.9, by=0.1)), iterations=100) {
  #Get current sample sizes
  size.sample <- colSums(count.table)
  #Get overall taxa count
  n.tax <- nrow(count.table)
  
  #Get relative counts
  count.table <- as.matrix(count.table)
  ct.rel <- as.matrix(t(t(count.table) / size.sample))
  
  #Preallocate results collector data.frame
  #=> add trivial rarefaction step at 0 counts
  results.rarefy <- data.frame(
    sample.name = colnames(count.table),
    N.seq = 0,
    N.obs = 0
  )
  
  #Iterate through rarefaction steps and perform rarefactions
  for (step in steps) {
    curr.sizes <- ceiling(size.sample * step)
    curr.N_obs <- mapply(function(rel.counts, n.rare, n.tax) {
      tmp.counts <- replicate(iterations, sample(1:n.tax, size = n.rare, prob = rel.counts, replace = T))
      mean(apply(tmp.counts, 2, function(x) {length(unique(x))}))
    }, rel.counts=as.data.frame(ct.rel), n.rare=curr.sizes, MoreArgs = list(n.tax=n.tax))
    #Append
    results.rarefy <- rbind(results.rarefy, data.frame(
      sample.name = colnames(count.table),
      N.seq = curr.sizes,
      N.obs = curr.N_obs
    ))
  }
  
  #Add N.obs at 100% counts
  results.rarefy <- rbind(results.rarefy, data.frame(
    sample.name = colnames(count.table),
    N.seq = colSums(count.table),
    N.obs = colSums(count.table > 0)
  ))
  
  #Return
  results.rarefy
}
################################################################################
################################################################################

