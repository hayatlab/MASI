de <- function(df,test='MAST'){
  library(Seurat)
  library(Matrix)
  celltype <- as.factor(df$celltype)
  df <- df[,1:(ncol(df)-1)]
  df <- t(df)
  df <- as(df, "dgTMatrix")
  colnames(df)<-paste('cell',c(1:ncol(df)),sep = "")
  names(celltype)<-paste('cell',c(1:ncol(df)),sep = "")
  pbmc <- CreateSeuratObject(counts = df, project = "source", min.cells = 0, min.features = 0)
  
  pbmc@active.ident <- celltype
  pbmc$seurat_clusters <- celltype
  
  diff_genes<-FindAllMarkers(pbmc,test.use=test,only.pos = TRUE,verbose = FALSE,min.cells.group = 2)
  return(data.frame(diff_genes))
}