agreement <- function(Beta, patientid) {
  pt.table <- table(patientid)
  sum.over.pt <- 0
  for (j in (1:length(pt.table)))  {
    if(pt.table[j]==1) {next()}
    which.pt <- which(patientid == names(pt.table)[j])
    mat.temp <- Beta[,which.pt]
    mat.diffs <- apply(mat.temp, 1,
                       function(x) { tmp <- outer(x,x,"-")^2
                       sum(tmp[upper.tri(tmp)])/
                         sum(upper.tri(tmp) )   } )
    
    sum.over.pt <- sum.over.pt + mat.diffs
  } # end for
  return(mean(sum.over.pt, na.rm=TRUE))
}
