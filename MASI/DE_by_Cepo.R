de <- function(df){
  library(Cepo)

  celltype <- as.factor(df$celltype)
  df <- df[,1:(ncol(df)-1)]
  df <- t(df)
  colnames(df)<-paste('cell',c(1:ncol(df)),sep = "")
  names(celltype)<-paste('cell',c(1:ncol(df)),sep = "")
  
  ds_res = Cepo(exprsMat = df,cellType = celltype)
  res_name = topGenes(object = ds_res, n = 20)
  
  return(data.frame(res_name))
}