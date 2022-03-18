rra <- function(df){
  library(RobustRankAggreg)
  glist <- list()
  for (i in c(1:ncol(df))) {
    glist[[i]]<-as.character(df[,i])
  }
  ranks<-aggregateRanks(glist = glist, N = 50, method = "stuart")
  return(data.frame(ranks))
}