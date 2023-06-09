# MASI

MASI: marker-assisted standardization and integration for single-cell transcriptomics data

Manuscript: Fast model-free standardization and integration of single-cell transcriptomics data

### 1. Brief description
MASI utilizes robust marker idenfication to identify marker genes from reference data and transfers cell-type labels to target data through MACA.

![alt text](https://github.com/hayatlab/MASI/blob/main/MASI/Figure%201.jpg?raw=true)

### 2. Install requirement packages
    pip install scanpy cosg rpy2 
    pip install sccaf ##optional
    pip install fa2 ##install if doing integrative lineage analysis
    
    ##install Seurat and RobustRankAggreg separately in R
    install.packages('Seurat')
    install.packages('RobustRankAggreg')
    
    ##Installment of PyTorch (Optional)
    ##We noticed adding BatchNorm1D as data transformation can further remove batch effects 
    ##but may sacrifice discriminative power for cell-type identification.
    
    ##Then git clone MASI to the working directory for use.
    
### 3. Usage
    import MASI as masi
    
    ##identify cell-type markers from reference data
    ##source is Anndata of reference data that has cell-type labels named 'cell_type'
    cell_markers = masi.marker_identification(source_data=source)
    
    ##step 1 transform gene expression matrix to cell-type score matrix
    ##ad is combined expression data in Anndata format, and cell_markers is dict for cell-type markers 
    ##scores can further be used for visualization and other downstream analyses
    scores, labels = masi.gene2cell(ad=ad,cell_markers=cell_markers,use_weight=True)
    ##step 2 clustering and parallel annotation
    annotation= masi.parallel(scores=scores,labels=labels,batch_size=50000)

### 4. Reproduce results in manuscript
Please see tutorials at https://github.com/hayatlab/MASI/tree/main/tutorial

Processed data can be found at https://figshare.com/articles/dataset/Fast_model-free_integration_and_transfer_learning_via_MASI_for_single-cell_expression_data/18866264.
    
### 5. Citation
Xu et al. "MASI enables fast model-free standardization and integration of single-cell transcriptomics data". <a href="https://www.nature.com/articles/s42003-023-04820-3">Communications Biology</a>

### 6. Integration of multimodal single-cell and spatial data
A major update of MASI is availale at https://github.com/hayatlab/MASIv2

MASIv2 integrates large-scale multimodal single-cell and spatial data with few computational resources. <a href="https://www.biorxiv.org/content/10.1101/2023.05.15.540808v1">bioRxiv</a>
