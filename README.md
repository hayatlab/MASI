# MASI

MASI: marker-assisted standardization and integration for single-cell transcriptomics data

Manuscript: Fast model-free standardization and integration of single-cell transcriptomics data

### 1. Brief description
MASI utilizes robust marker idenfication to identify marker genes from reference data and transfers cell-type labels to target data through MACA.

![alt text](https://github.com/hayatlab/MASI/blob/main/MASI/Figure%201.jpg?raw=true)

### 2. Install requirement packages
    pip install scanpy cosg rpy2 sccaf
    pip install fa2##install if doing integrative lineage analysis
    
    ##install Seurat and RobustRankAggreg separately in R
    install.packages('Seurat')
    install.packages('RobustRankAggreg')
    
    ##Installment of PyTorch (Optional)
    ##We noticed adding BatchNorm1D as data transformation can further remove batch effects 
    ##but may sacrifice discriminative power for cell-type identification.
    
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
    

### 5. References
    1. Xu, Y., et al. "MACA: Marker-based automatic cell-type annotation for single cell expression data." Bioinformatics (2021).
    2. Wolf, F. Alexander, Philipp Angerer, and Fabian J. Theis. "SCANPY: large-scale single-cell gene expression data analysis." Genome biology (2018).
    3. Butler, Andrew, et al. "Integrating single-cell transcriptomic data across different conditions, technologies, and species." Nature biotechnology (2018).

### 6. Citation
Xu et al. "Fast model-free standardization and integration of single-cell transcriptomics data". <a href="https://www.biorxiv.org/content/10.1101/2022.03.28.486110v1">bioRxiv</a>

## update 05/18/2022
We are working on a major update of MASI to integrate single-cell multi-omics data.
