# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:31:17 2020

@author: Yang Xu
"""

import umap
import scipy
import anndata
import collections
import numpy as np
import pandas as pd
import cosg as cosg
import scanpy as sc
import multiprocessing

#import torch
#import torch.nn.functional as F

import rpy2.robjects as ro
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from sklearn.metrics import confusion_matrix
from sklearn.neighbors import KNeighborsClassifier
from sklearn.feature_extraction.text import TfidfTransformer

from SCCAF import SCCAF_optimize_all

import warnings
warnings.filterwarnings("ignore")

##-----------------------------------------------------------------------------
##main functions of transfer learning via parallel computing
def ensemble_labels(multi_labels=None):
    ensemble = []
    for i in range(multi_labels.shape[0]):
        ks=[]
        vs=[]
        for k,v in collections.Counter(multi_labels[i,:].tolist()).items():
            ks.append(k)
            vs.append(v)
        ensemble.append(ks[vs.index(max(vs))])
    return ensemble

#def BatchNorm(ad):
    ##use torch batchnorm1d layer to perform batch normalization
    ##not used in final method
#    X = ad.X.copy()
#    if scipy.sparse.issparse(X):
#        X = X.todense()
    
#    X_all_tensor = torch.tensor(X).float()
#    batchnorm = torch.nn.BatchNorm1d(X.shape[1], affine=False)
#    activation = torch.nn.LeakyReLU(0.2)
#    batchnorm.to(torch.device("cpu"))
#    activation.to(torch.device("cpu"))

#    newX = np.zeros((X.shape))

    #batches = list(set(ad.obs['study'].values.tolist()))
    #for b in batches:
    #    pred = batchnorm(X_all_tensor[ad.obs['study']==b,:])
    #    pred = activation(pred)
    #    #pred = F.normalize(pred,p=1,dim=1)
    #    pred = torch.Tensor.cpu(pred).detach().numpy()
    #    #scaler = StandardScaler()
    #    #pred=scaler.fit_transform(X[ad.obs['study']==b,:])
    #    newX[ad.obs['study']==b,:]=pred

#    batch_size=1024
#    for j in range(X.shape[0]//batch_size+1):
#        pred = batchnorm(X_all_tensor[j*batch_size:(j+1)*batch_size,:])
#        pred = activation(pred)
#        pred = F.normalize(pred, dim=1,p=2)
#        pred = torch.Tensor.cpu(pred).detach().numpy()
#        newX[j*batch_size:(j+1)*batch_size,:]=pred
    
#    ad.obsm['X_batch']=newX

#    return ad

def gene2cell(ad=None, cell_markers=None,use_weight=False,thresh=0.25,
              if_tfidf=True,if_thresh=True,use_knn=False,use_umap=True):
    ##TF-IDF transformation
    X = ad.X.copy()
    if scipy.sparse.issparse(X):
        X = X.todense()
    
    if if_tfidf:
        tf_transformer = TfidfTransformer(use_idf=True).fit(X)
        X= tf_transformer.transform(X).todense()
    
    labels = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    labels.columns = cell_markers.keys()
    exprsed = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    exprsed.columns = cell_markers.keys()
    celltype_size = {}
    
    ##create artifical labels for each cell
    if use_weight == True:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            marker_index = -1
            for i in v:
                marker_index += 1
                if i in ad.var.index:
                    if if_thresh:
                        expr95 = np.percentile(X[:,ad.var.index == i],95)
                        thresh = thresh * expr95
                        l = np.array(X[:,ad.var.index == i])
                        l[X[:,ad.var.index == i]<=thresh]=0
                    else:
                        l = np.array(X[:,ad.var.index == i])
                    ##consider marker weight
                    l = l*(1-marker_index/(len(v)*2))##default 2
                    
                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0]) 

    else:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            for i in v:
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = thresh * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0
                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0])        
    
    if use_knn:##not used in final method
        if use_umap:
            down_samp = pd.DataFrame(labels.values[ad.obs['source']=='reference',:])
            umaps = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2,
                              metric="cosine").fit(down_samp.values)
            embedding = umaps.transform(labels.values)
            ad.obsm['X_score']=embedding
        else:
            ad.obsm['X_score']=labels.values
        subsample = ad[ad.obs['source']=='reference']
        subsample = downsample_to_smallest_category(subsample, 'cell_type', min_cells=1000, 
                                                    keep_small_categories=True)#default 500
        neigh = KNeighborsClassifier(n_neighbors=5,weights='distance')
        neigh.fit(subsample.obsm['X_score'], subsample.obs['cell_type'])
        new_labels = neigh.predict_proba(ad.obsm['X_score'])
        new_labels = pd.DataFrame(new_labels)
        new_labels.columns = neigh.classes_
        
        labels = labels.reindex(sorted(labels.columns), axis=1)
        new_labels = new_labels.reindex(sorted(new_labels.columns), axis=1)
    
    else:
        assess1 = np.argmax((labels*exprsed).values,axis=1)
        vals1 = 0
        for k,v in collections.Counter(assess1).items():
            if v >= 5:
                vals1 += 1
                            
        assess1 = vals1
    
        assess2 = np.argmax((labels).values,axis=1)
        vals2 = 0
        for k,v in collections.Counter(assess2).items():
            if v >= 5:
                vals2 += 1
                       
        assess2 = vals2
    
        assess = [assess1,assess2]
    
        new_labels = [labels*exprsed,labels][assess.index(max(assess))]
        #new_labels = labels*exprsed##consider the number of expressed marker of each cell-type for each cell
        #new_labels = labels#*exprsed
    
    print(labels.shape)
    return labels, new_labels

def multiMASI(labels=None,new_labels=None):
    
    ad = anndata.AnnData(X=labels)
    ##create Label1
    labels1 = np.argmax(new_labels.values,axis=1)
    ad.obs['Label1'] = labels1
    
    ad.obsm['Score']=labels.values
    
    ##two key parameters for louvain clustering
    res = [3,5,7]
    n_neis = [5,10,15]
    
    label_list = np.zeros((ad.X.shape[0],len(res)*len(n_neis))).astype('str')
    indexs = 0
    for r in res:
        for nei in n_neis:
            
            sc.pp.neighbors(ad, use_rep="Score", n_neighbors=nei,metric='cosine')
            sc.tl.louvain(ad, resolution=r, key_added = 'louvain')
            #sc.tl.leiden(ad, resolution=r, key_added = 'louvain')
            
            cm = confusion_matrix(ad.obs['louvain'].values.astype(int),
                                  labels1)
            
            normed_cm = cm.copy().T
            normed_cm = normed_cm/np.sum(normed_cm,axis=0)
            normed_cm = np.nan_to_num(normed_cm)
            normed_cm = normed_cm.T
            mapping={}
            
            mapping = np.argmax(normed_cm,axis=1)
            mapmax = np.max(normed_cm,axis=1)
            mapmax = np.nan_to_num(mapmax)
            
            clustering = ad.obs['louvain'].values.astype(int)
            new_cluster = np.zeros((len(clustering)))
            for i in range(len(mapping)):
                tof = clustering==i
                if mapmax[i]>=0.5:
                    new_cluster[tof]=mapping[i]
                else:
                    sub_labels = new_labels.values[tof,:]
                    if sub_labels.shape[0]>0:
                        freqs = []
                        for j in range(sub_labels.shape[0]):
                            zscore=scipy.stats.zscore(sub_labels[j,:])
                            orderi = np.array([g for g in range(sub_labels.shape[1])])[zscore>3].tolist()
                            orderj = np.argsort(sub_labels[j,:])[::-1][:3].tolist()#default 3
                            a = [len(orderi),len(orderj)]
                            freqs += [orderi,orderj][a.index(max(a))]
                        vals = 0
                        ks = 0
                        for k,v in collections.Counter(freqs).items():
                            if v >= vals:
                                ks = k
                                vals = v
                        
                        if vals/sub_labels.shape[0]>=0.5:
                            new_cluster[tof]=ks
                        else:
                            new_cluster[tof]=-1
    
            mapped =new_cluster.astype('int')
            mapped=mapped.astype('str')
            ad.obs['Mapped'] = mapped
            
            cell_dict = {}
            for k,v in collections.Counter(new_cluster.tolist()).items():
                if int(k)>=0:
                    cell_dict[int(k)]=new_labels.columns[int(k)]
                else:
                    cell_dict[int(k)]="unassigned"
            label_list[:,indexs]=ad.obs['Mapped'].values
            indexs+=1
            
    return label_list

def parallel(scores=None,labels=None,batch_size=20000,n_core=10):
    index = np.array([i for i in range(scores.shape[0])])
    
    r = np.random.permutation(scores.shape[0])
    r_index = index[r]
    r_scores = scores.iloc[r,:]
    r_label1 = labels.iloc[r,:]
    scores_list= []
    for j in range(scores.shape[0]//batch_size+1):
        scores_list.append((r_scores.iloc[j*batch_size:(j+1)*batch_size,:],
                            r_label1.iloc[j*batch_size:(j+1)*batch_size,:]))
            
    pool = multiprocessing.Pool(processes=n_core)
    mapped = pool.starmap(multiMASI, scores_list)

    merged = mapped[0]
    for m in mapped[1:]:
        merged = np.concatenate((merged,m),0)
    merged = merged[r_index.argsort(),:]
    
    ensemble = ensemble_labels(merged)
    ensemble = np.array(ensemble)
    
    annotations=[]
    for e in ensemble:
        if int(e)>=0:
            annotations.append(labels.columns[int(e)])
        else:
            annotations.append("unassigned")
        
    return annotations

##-----------------------------------------------------------------------------
##robust marker identification
def downsample_to_smallest_category(adata,column="cell_type",random_state=None,
                                    min_cells=15,keep_small_categories=False):
    
    counts = adata.obs[column].value_counts(sort=False)
    min_size = min(counts[counts >= min_cells])
    sample_selection = None
    for sample, num_cells in counts.items():
        if num_cells <= min_cells:
            if keep_small_categories:
                sel = adata.obs.index.isin(
                    adata.obs[adata.obs[column] == sample].index)
            else:
                continue
        else:
            sel = adata.obs.index.isin(
                adata.obs[adata.obs[column] == sample]
                .sample(min_size, random_state=random_state)
                .index
            )
        if sample_selection is None:
            sample_selection = sel
        else:
            sample_selection |= sel
    return adata[sample_selection].copy()

def marker_identification_fast(source_data=None,diff_method="wilcoxon",
                               repeats=10,num_sub_cells=500,if_metacells=False):
    all_cell_types = list(set(source_data.obs['cell_type']))
    cell_markers = []
    for i in range(repeats):#default 50
        if if_metacells:
            subsample = source_data.copy()
            subsample = downsample_to_smallest_category(subsample, 'cell_type', min_cells=50, 
                                                        keep_small_categories=True)
        else:
            subsample = source_data.copy()
            subsample = downsample_to_smallest_category(subsample, 'cell_type', min_cells=num_sub_cells, 
                                                        keep_small_categories=True)#default 500
        if diff_method == 'cosg':
            cosg.cosg(subsample,key_added='cosg',mu=1,n_genes_user=50,groupby='cell_type')#,expressed_pct=0.5
            cellmarkers = pd.DataFrame(subsample.uns['cosg']['names']).iloc[:20,:]#default 50
        elif diff_method in ['t-test','t-test_overestim_var','wilcoxon','logreg']:
            sc.tl.rank_genes_groups(subsample, 'cell_type', method=diff_method)#,tie_correct=True)##Seurat wilcoxon
            cellmarkers = pd.DataFrame(subsample.uns['rank_genes_groups']['names']).iloc[:20,:]#default 50
        else:
            r = robjects.r
            r['source']('DE_by_Seurat.R')
            de_r = robjects.globalenv['de']
            
            X = subsample.X.copy()
            if scipy.sparse.issparse(X):
                X = X.todense()

            exprs = pd.DataFrame(X)
            exprs.columns = subsample.var.index.tolist()
            exprs['celltype']=subsample.obs['cell_type'].values
            
            with localconverter(ro.default_converter + pandas2ri.converter):
                exprs_r = ro.conversion.py2rpy(exprs)
                
            result_r = de_r(exprs_r,test=diff_method)
            with localconverter(ro.default_converter + pandas2ri.converter):
                result = ro.conversion.rpy2py(result_r)

            cellmarkers={}
            allcelltypes = list(set(subsample.obs['cell_type']))
            for f in allcelltypes:
                celltypes = result[result['cluster']==f]
                celltypes = celltypes.iloc[:20,:]['gene'].values.tolist()#default 20
                cellmarkers[f]=celltypes
            #cellmarkers = pd.DataFrame(cellmarkers)
            cellmarkers = pd.DataFrame.from_dict(cellmarkers, orient='index').T
            cellmarkers = cellmarkers.fillna('')
        cell_markers.append(cellmarkers)
        
    cell_markers_db = {}
    for c in all_cell_types:
        cell_markers_db[c]=[]
        for i in cell_markers:
            cell_markers_db[c].append(i[c].values.tolist())
    
    ##r script of Robust rank aggregation
    r = robjects.r
    r['source']('RobustRankAggreg.R')
    rra_r = robjects.globalenv['rra']
    
    cell_markers_score ={}
    cell_markers_rank ={}
    for c in all_cell_types:
        celltype = pd.DataFrame(cell_markers_db[c]).T
        with localconverter(ro.default_converter + pandas2ri.converter):
            celltype_r = ro.conversion.py2rpy(celltype)
        
        result_r = rra_r(celltype_r)
        with localconverter(ro.default_converter + pandas2ri.converter):
            result = ro.conversion.rpy2py(result_r)
            
        result = result[result['Score'].values<=0.9].iloc[:20,:]
        cell_markers_rank[c]=result['Name'].values.tolist()
        cell_markers_score[c]=result['Score'].values.tolist()
    
    return cell_markers_rank, cell_markers_score

def marker_identification(source_data=None,num_sub_cells=500,if_metacells=False):
    
    all_cell_types = list(set(source_data.obs['cell_type']))
    cell_markers = []
    
    if if_metacells:
        #repeats = ['cosg','t-test_overestim_var','wilcoxon-tie','bimod','wilcox','roc','MAST']
        #repeats = ['MAST','bimod','roc','t-test_overestim_var','wilcoxon-tie','cosg']
        repeats = ['poisson','bimod','roc','cosg','t-test_overestim_var','wilcoxon-tie']
    else:
        #repeats = ['poisson','cosg','bimod']
        repeats = ['poisson','bimod','roc','MAST','wilcox','cosg','t-test_overestim_var','wilcoxon-tie']#
    #repeats = ['poisson','cosg','bimod']##option 1
    #repeats = ['cosg','t-test_overestim_var','bimod','poisson']##option 2 and it works best so far
    #repeats = ['cosg','t-test_overestim_var','bimod','poisson','wilcox','negbinom','MAST']##option 3
    for i in repeats:
        if if_metacells:
            subsample = source_data.copy()
            subsample = downsample_to_smallest_category(subsample, 'cell_type', min_cells=50, 
                                                        keep_small_categories=True)
        else:
            subsample = source_data.copy()
            subsample = downsample_to_smallest_category(subsample, 'cell_type', min_cells=num_sub_cells, 
                                                        keep_small_categories=True)#default 500
        if i == 'cosg':
            cosg.cosg(subsample,key_added='cosg',mu=1,n_genes_user=50,groupby='cell_type')#,expressed_pct=0.5
            cellmarkers = pd.DataFrame(subsample.uns['cosg']['names']).iloc[:20,:]#default 50
        elif i == 'wilcoxon-tie':
            sc.tl.rank_genes_groups(subsample, 'cell_type', method='wilcoxon',tie_correct=True)##Seurat wilcoxon
            cellmarkers = pd.DataFrame(subsample.uns['rank_genes_groups']['names']).iloc[:20,:]#default 50
            print('**finished identifying marker genes by %s**' % i)
        elif i in ['t-test','t-test_overestim_var','logreg','wilcoxon']:
            sc.tl.rank_genes_groups(subsample, 'cell_type', method=i)#,tie_correct=True)##Seurat wilcoxon
            cellmarkers = pd.DataFrame(subsample.uns['rank_genes_groups']['names']).iloc[:20,:]#default 50
            print('**finished identifying marker genes by %s**' % i)
        else:
            #Seurat r script
            r = robjects.r
            r['source']('DE_by_Seurat.R')
            de_r = robjects.globalenv['de']
            
            X = subsample.X.copy()
            if scipy.sparse.issparse(X):
                X = X.todense()

            exprs = pd.DataFrame(X)
            exprs.columns = subsample.var.index.tolist()
            exprs['celltype']=subsample.obs['cell_type'].values
            
            with localconverter(ro.default_converter + pandas2ri.converter):
                exprs_r = ro.conversion.py2rpy(exprs)
                
            result_r = de_r(exprs_r,test=i)
            with localconverter(ro.default_converter + pandas2ri.converter):
                result = ro.conversion.rpy2py(result_r)

            cellmarkers={}
            allcelltypes = list(set(subsample.obs['cell_type']))
            for f in allcelltypes:
                celltypes = result[result['cluster']==f]
                celltypes = celltypes.iloc[:20,:]['gene'].values.tolist()#default 20
                cellmarkers[f]=celltypes
            #cellmarkers = pd.DataFrame(cellmarkers)
            cellmarkers = pd.DataFrame.from_dict(cellmarkers, orient='index').T
            cellmarkers = cellmarkers.fillna(' ')
            print('**finished identifying marker genes by %s**' % i)
        cell_markers.append(cellmarkers)
    
    cell_markers_db = {}
    for c in all_cell_types:
        cell_markers_db[c]=[]
        for i in cell_markers:
            cell_markers_db[c].append(i[c].values.tolist())
    
    ##r script of Robust rank aggregation
    r = robjects.r
    r['source']('RobustRankAggreg.R')
    rra_r = robjects.globalenv['rra']
    
    #cell_markers_score ={}
    cell_markers_rank ={}
    for c in all_cell_types:
        celltype = pd.DataFrame(cell_markers_db[c]).T
        with localconverter(ro.default_converter + pandas2ri.converter):
            celltype_r = ro.conversion.py2rpy(celltype)
        
        result_r = rra_r(celltype_r)
        with localconverter(ro.default_converter + pandas2ri.converter):
            result = ro.conversion.rpy2py(result_r)
            
        genes = result['Name'].values.tolist()
        newgenes = [g for g in genes if g != ' '][:20]
        result = result.iloc[:20,:]#[result['Score'].values<=0.9]
        cell_markers_rank[c]=newgenes
        #cell_markers_score[c]=result['Score'].values.tolist()
    
    return cell_markers_rank#, cell_markers_score

def marker_identification_LC(source_data=None,num_sub_cells=500,if_metacells=False):
    
    all_cell_types = list(set(source_data.obs['cell_type']))
    cell_markers = []
    
    if if_metacells:
        #repeats = ['cosg','t-test_overestim_var','wilcoxon-tie','bimod','wilcox','roc','MAST']
        #repeats = ['MAST','bimod','roc','t-test_overestim_var','wilcoxon-tie','cosg']
        repeats = ['MAST','poisson','bimod','wilcox','t-test_overestim_var','wilcoxon-tie']
    else:
        repeats = ['poisson','bimod','wilcox','t-test_overestim_var','wilcoxon-tie']#
    
    for i in repeats:
        if if_metacells:
            subsample = source_data.copy()
            subsample = downsample_to_smallest_category(subsample, 'cell_type', min_cells=50, 
                                                        keep_small_categories=True)
        else:
            subsample = source_data.copy()
            subsample = downsample_to_smallest_category(subsample, 'cell_type', min_cells=num_sub_cells, 
                                                        keep_small_categories=True)
        if i == 'wilcoxon-tie':
            sc.tl.rank_genes_groups(subsample, 'cell_type', method='wilcoxon',tie_correct=True)##Seurat wilcoxon
            cellmarker = pd.DataFrame(subsample.uns['rank_genes_groups']['names']).iloc[:50,:]
            cellpval = pd.DataFrame(subsample.uns['rank_genes_groups']['pvals']).iloc[:50,:]
            
            cellmarkers={}
            allcelltypes = list(set(subsample.obs['cell_type']))
            for f in allcelltypes:
                genes = cellmarker[f].values.tolist()
                pvals = cellpval[f].values.tolist()
                cellmarkers[f]=pd.DataFrame({'genes':genes,'Pval':pvals})
            print('**finished identifying marker genes by %s**' % i)
            
        elif i in ['t-test','t-test_overestim_var','logreg','wilcoxon']:
            sc.tl.rank_genes_groups(subsample, 'cell_type', method=i)#,tie_correct=True)##Seurat wilcoxon
            cellmarker = pd.DataFrame(subsample.uns['rank_genes_groups']['names']).iloc[:50,:]
            cellpval = pd.DataFrame(subsample.uns['rank_genes_groups']['pvals']).iloc[:50,:]
            
            cellmarkers={}
            allcelltypes = list(set(subsample.obs['cell_type']))
            for f in allcelltypes:
                genes = cellmarker[f].values.tolist()
                pvals = cellpval[f].values.tolist()
                cellmarkers[f]=pd.DataFrame({'genes':genes,'Pval':pvals})
            print('**finished identifying marker genes by %s**' % i)
            
        else:
            #Seurat r script
            r = robjects.r
            r['source']('DE_by_Seurat.R')
            de_r = robjects.globalenv['de']
            
            X = subsample.X.copy()
            if scipy.sparse.issparse(X):
                X = X.todense()

            exprs = pd.DataFrame(X)
            exprs.columns = subsample.var.index.tolist()
            exprs['celltype']=subsample.obs['cell_type'].values
            
            with localconverter(ro.default_converter + pandas2ri.converter):
                exprs_r = ro.conversion.py2rpy(exprs)
                
            result_r = de_r(exprs_r,test=i)
            with localconverter(ro.default_converter + pandas2ri.converter):
                result = ro.conversion.rpy2py(result_r)

            cellmarkers={}
            allcelltypes = list(set(subsample.obs['cell_type']))
            for f in allcelltypes:
                celltypes = result[result['cluster']==f]
                celltype = celltypes['gene'].values.tolist()[:50]#default 20
                cellpval = celltypes['p_val'].values.tolist()[:50]#default 20
                cellmarkers[f]=pd.DataFrame({'genes':celltype,'Pval':cellpval})
                
            print('**finished identifying marker genes by %s**' % i)
        cell_markers.append(cellmarkers)
    
    ##r script of lancaster combination
    r = robjects.r
    r['source']('LancasterCombination.R')
    LC = robjects.globalenv['lancaster.combination']
    
    #cell_markers_score ={}
    cell_markers_rank ={}
    for c in all_cell_types:
        celltype = cell_markers[0][c]
        celltype.index = celltype['genes'].values
        for de in cell_markers[1:]:
            celltype2 = de[c]
            celltype2.index = celltype2['genes'].values
            celltype = celltype.merge(celltype2, on='genes', how='inner')
        celltype.index = celltype['genes'].values
        celltype.pop('genes')
            
        with localconverter(ro.default_converter + pandas2ri.converter):
            celltype_r = ro.conversion.py2rpy(celltype)
        
        result_r = LC(celltype_r)
        with localconverter(ro.default_converter + pandas2ri.converter):
            result = ro.conversion.rpy2py(result_r)
            
        genes = result.index.tolist()[:20]
        cell_markers_rank[c]=genes
        #cell_markers_score[c]=result['Score'].values.tolist()
    
    return cell_markers_rank#, cell_markers_score

##-----------------------------------------------------------------------------
##metacells via louvain community detection
def metacells(source_data=None,res=10):
    
    sc.tl.pca(source_data, svd_solver='arpack')
    sc.pp.neighbors(source_data, n_neighbors=5,metric='cosine')
    sc.tl.louvain(source_data, resolution=res, key_added = 'louvain')
    
    
    obs_codes, obs_names = pd.factorize(source_data.obs['cell_type'])
    pred_codes, pred_names = pd.factorize(source_data.obs['louvain'])
    cm = confusion_matrix(obs_codes, pred_codes)
    norm_cm = (cm.T / cm.sum(axis=1)).T
    norm_cm = pd.DataFrame(norm_cm)
    norm_cm = norm_cm.iloc[:len(obs_names),:len(pred_names)]
    norm_cm.index = obs_names
    norm_cm.columns = pred_names
    norm_cm = norm_cm[sorted(norm_cm.columns)]
    norm_cm = norm_cm.sort_index(axis=0)
    
    louvain2celltype = np.argmax(norm_cm.values, axis=0)
    meta = {}
    for i in range(norm_cm.shape[1]):
        meta[norm_cm.columns[i]]=norm_cm.index.tolist()[louvain2celltype[i]]
    meta = pd.DataFrame.from_dict(meta,orient='index')
    meta.columns = ['cell_type']
    
    unique_keys, row = np.unique(source_data.obs['louvain'], return_inverse=True)
    x ={}
    if scipy.sparse.issparse(source_data.X):
        for k in unique_keys:
            x[k]=source_data[source_data.obs['louvain']==k].X.mean(axis=0).reshape(source_data.X.shape[1])[0].tolist()[0]
    else:
        for k in unique_keys:
            x[k]=np.mean(source_data[source_data.obs['louvain']==k].X,axis=0)
    x = pd.DataFrame.from_dict(x,orient='index')
    x.columns = source_data.var.index
    metacells = anndata.AnnData(X=x,obs=meta)
    
    return metacells

##-----------------------------------------------------------------------------
##Use SCCAF to identify subtypes, when reference data has less annotation 
##resolution than target data (Miao et al., Nature Methods, 2020)
def subtype_identification(ad=None,cell_type_score=None):
    
    ad.obsm['X_pca']=cell_type_score.values
    sc.pp.neighbors(ad, use_rep="X_pca", n_neighbors=5,metric='cosine')
    
    all_anns = list(set(ad.obs['Annotation'].values.tolist()))
    ad.obs['newAnnotation']=" "
    newAnnotation = np.array(ad.obs['newAnnotation'])
    for a in all_anns:
        ad_sccaf = ad[ad.obs['Annotation']==a]
        sc.tl.leiden(ad_sccaf, resolution=0.2, key_added='L2_Round0')
        SCCAF_optimize_all(min_acc=0.9, ad=ad_sccaf, basis ='umap', use='pca',
                           prefix = 'L2')
        subcluster = np.array(ad_sccaf.obs['L2_result'].values.tolist())
        subid = pd.DataFrame.from_dict(collections.Counter(ad_sccaf.obs['L2_result'].values.tolist()),
                                       orient='index')
        subid = subid[subid[0]<=15]
        for i in subid.index.tolist():
            subcluster[subcluster==i]=0
        newAnnotation[ad.obs['Annotation']==a]=subcluster#ad_sccaf.obs['L2_result'].values.tolist()
        
    subtypes = []
    for i in range(len(newAnnotation)):
        subtypes.append(ad.obs['Annotation'].values[i]+"-"+newAnnotation[i])
    ad.obs['newAnnotation']=subtypes
    
    return ad