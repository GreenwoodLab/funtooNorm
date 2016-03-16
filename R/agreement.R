# Function to measure intra-replicate agreement in methylation data.
agreement <- function(Beta, individualID) {
  pt.table <- table(individualID)
  sum.over.pt <- 0
  for (j in (1:length(pt.table)))  {
    if(pt.table[j]==1) {next()}
    which.pt <- which(individualID == names(pt.table)[j])
    mat.temp <- Beta[,which.pt]
    combs=combn(1:ncol(mat.temp),2)
    mat.diffs=rowSums(apply(combs,2,function(x){(mat.temp[,x[1]]-mat.temp[,x[2]])^2}))/ncol(combs)
    
    sum.over.pt <- sum.over.pt + mat.diffs
  } 
  return(mean(sum.over.pt, na.rm=TRUE))
}
