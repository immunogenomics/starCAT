#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import scanpy as sc
import os
import scipy.sparse as sp
import yaml
import requests
import tarfile
import warnings
# from starcat import load_df_from_npz

def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj


def df_col_corr(X, Y):
    '''
    Compute pairwise Pearson correlation matrix of columns of X and Y. Returns
    R which is n_X_cols X n_Y_cols
    '''

    if X.isnull().any().any() or Y.isnull().any().any():
        raise Exception("NaNs found in spectra matrix which efficient pearson correlation does not support")

    X_norm = X.subtract(X.mean(axis=0), axis=1)
    X_norm= X_norm.divide(X_norm.std(axis=0), axis=1)
    Y_norm = Y.subtract(Y.mean(axis=0), axis=1)
    Y_norm= Y_norm.divide(Y_norm.std(axis=0), axis=1)
    R = np.dot(X_norm.T, Y_norm) / (X.shape[0]-1)
    
    R = pd.DataFrame(R, index=X.columns, columns=Y.columns)
    return(R)

cnmf_dir_strs = {
                'nmf_genes_list' :  os.path.join('{odir}', '{name}', '{name}'+'.overdispersed_genes.txt'),
                'tpm_stats' :  os.path.join('{odir}', '{name}', 'cnmf_tmp', '{name}'+'.tpm_stats.df.npz'),
                'gene_spectra_score__txt': os.path.join('{odir}', '{name}', '{name}'+'.gene_spectra_score.k_{k}.dt_{ldstr}.txt'),
                'gene_spectra_tpm__txt': os.path.join('{odir}', '{name}', '{name}'+'.gene_spectra_tpm.k_{k}.dt_{ldstr}.txt'),
                'starcat_spectra__txt': os.path.join('{odir}', '{name}', '{name}'+'.starcat_spectra.k_{k}.dt_{ldstr}.txt'),
            }


class BuildConsensusReference():  
    def __init__(self, cnmf_paths, output_dir='.', prefix = '', ks = None, density_thresholds = None,
                 tpm_fns = None, score_fns = None, order_thresh = None, corr_thresh = 0.5, pct_thresh = 0.666):
        """
        Class for building consensus gene expression programs (GEPs) from 2 or more cNMF results.
        
        Parameters
        ----------
        cnmf_paths : list, paths to cNMF project directories for clustering. Should include the cNMF project
            name at the end of the path, I.e. [cnmf_output_dir]/[cnmf_name]

        output_dir : path, optional (default=".") output directory for spectra results

        prefix: str, prefix to add to the output filenames (default='')

        ks: list, list of K (integers) used for cNMF results

        density_thresholds: list, list of density thresholds (floats) used for cNMF results

        tpm_fns : list, optional list of paths to cNMF TPM spectra paths to cluster (use instead of ks, dts)

        score_fns : list, optional list of paths to cNMF spectra score paths (use instead of ks, dts)
        
        order_thresh: int, maximum rank for 2 programs to be in each other's top correlated programs and allowed
            to cluster (default=len(cnmf_objs))

        corr_thresh: float, minimum correlation for programs to cluster
        
        pct_thresh: float, minimum number of connected programs to add a program to a cluster or merge clusters 

        """
        
        # Load filepaths from cnmf object
        cnmf_paths = [os.path.normpath(x) for x in cnmf_paths]
        odir_paths = []
        self.dataset_names = []
        for path in cnmf_paths:
            if not os.path.exists(path):
                raise Exception("Input path %s does not exist" % path)

            npath = os.path.dirname(path)
            name = os.path.basename(path)
            odir_paths.append(npath)
            self.dataset_names.append(name)
            
        if (any([ks, density_thresholds])) & (not any([tpm_fns, score_fns])):
            tpm_fns = []
            score_fns = []
            for i in range(len(cnmf_paths)):
                ldstr= str(density_thresholds[i]).replace('.', '_')
                tpm_fns.append(cnmf_dir_strs['gene_spectra_tpm__txt'].format(odir=odir_paths[i], name=self.dataset_names[i], k=ks[i], ldstr=ldstr))
                score_fns.append(cnmf_dir_strs['gene_spectra_score__txt'].format(odir=odir_paths[i], name=self.dataset_names[i], k=ks[i], ldstr=ldstr))

        elif not (not any([ks, density_thresholds])) & (any([tpm_fns, score_fns])):
            raise TypeError('Object types %s, %s and %s, %s are not valid. Please pass only ks/density_thresholds OR \
            tpm_fns/score_fns.' % (type(ks), type(density_thresholds), type(tpm_fns), type(score_fns)))
        
        for i in range(len(cnmf_paths)):
            if any([cnmf_paths[i] != os.path.dirname(score_fns[i])]) or any([cnmf_paths[i] != os.path.dirname(tpm_fns[i])]):
                warnings.warn('At least one cnmf directory path is inconsistent with a spectra path.\
                              IF inputting score and tpm files, ensure they are in the same order as cnmf_paths')
        
        self.num_results = len(cnmf_paths)    
        self.spectra_tpm_all = []
        self.spectra_score_all = []
        self.hvgs_all = []
        self.stds_all = []
        
        self.output_dir = output_dir 
        self.prefix = prefix
        self.corr_thresh = corr_thresh
        self.pct_thresh = pct_thresh
        if order_thresh==None:
            self.order_thresh = self.num_results
        else:
            self.order_thresh = order_thresh

        # Load spectra, stats, HVGs
        for i, (tpm_fn, score_fn) in enumerate(zip(tpm_fns, score_fns)):
            tpm_stats = load_df_from_npz(cnmf_dir_strs['tpm_stats'].format(odir=odir_paths[i], name=self.dataset_names[i]))
            stds = tpm_stats['__std']
            spectra_tpm = pd.read_csv(tpm_fn, index_col = 0, sep = '\t')
            spectra_tpm = spectra_tpm.loc[:, ~spectra_tpm.columns.str.contains('AB_|prot')]
            spectra_tpm.index = '%s:' % (self.dataset_names[i]) + spectra_tpm.index.astype(str)
            spectra_score = pd.read_csv(score_fn, index_col = 0, sep = '\t')
            spectra_score.index = '%s:' % (self.dataset_names[i]) + spectra_score.index.astype(str)
            hvgs = open(cnmf_dir_strs['nmf_genes_list'].format(odir=odir_paths[i], name=self.dataset_names[i])).read().split('\n')
            
            self.spectra_tpm_all += [spectra_tpm]
            self.spectra_score_all += [spectra_score]
            self.hvgs_all += [hvgs]
            self.stds_all += [stds]
            

    def cluster_cnmf_results(self):
        """
        Perform and save clustering of cNMF results. 

        Parameters
        ----------

        """

        all_genes = list(set(spectra_tpm.columns) for spectra_tpm in self.spectra_tpm_all)
        intersect_genes_all = sorted(set.intersection(*all_genes))
        union_genes_all = sorted(set.intersection(*all_genes))
        
        # Renormalize and variance-normalize TPM spectra for all cNMF objects
        merged_data = {'TPM_Renorm_VarNorm':[], 'Scores':[] }

        for n in range(self.num_results):
            spectra_scores = self.spectra_score_all[n]
            spectra_renorm = self.spectra_tpm_all[n][intersect_genes_all].copy()
            spectra_renorm = spectra_renorm.div(spectra_renorm.sum(axis=1), axis=0)*1e6
            spectra_varnorm = spectra_renorm.div(self.stds_all[n][spectra_renorm.columns])
            new_genes = sorted(set(union_genes_all) - set(spectra_scores.columns))
            spectra_scores[new_genes] = np.nan
            spectra_scores = spectra_scores[union_genes_all]

            merged_data['TPM_Renorm_VarNorm'].append(spectra_varnorm)
            merged_data['Scores'].append(spectra_scores)

        merged_data['TPM_Renorm_VarNorm'] = pd.concat([x[intersect_genes_all] for x in merged_data['TPM_Renorm_VarNorm']], axis=0)        
        merged_data['Scores'] = pd.concat(merged_data['Scores'], axis=0)      

        # Define pairwise correlations for all GEPs
        self.correlate_geps()

        # Cluster GEPs using correlation matrix
        clus_df, clus_list = self.cluster_geps()

        spectra_tpm_grouped = merged_data['TPM_Renorm_VarNorm'].groupby(clus_list).mean()
        spectra_tpm_grouped.index = 'cGEP' + spectra_tpm_grouped.index.astype(str)
        spectra_scores_grouped = merged_data['Scores'].groupby(clus_list).mean()
        spectra_scores_grouped.index = 'cGEP' + spectra_scores_grouped.index.astype(str)
        hvgs_union = sorted(set.union(*list(set(hvgs) for hvgs in self.hvgs_all)).intersection(spectra_tpm_grouped))
        spectra_tpm_grouped_hvg = spectra_tpm_grouped[hvgs_union]

        # Filter singletons and save results
        filtsingle_index = clus_df[np.sum(-clus_df.isna(), axis = 1)>1].index
        clus_df.to_csv(os.path.join(self.output_dir, 
                                    '%sstarcat_consensus_clustering.txt' % self.prefix), '\t')
        clus_df.loc[filtsingle_index, :].to_csv(os.path.join(self.output_dir, 
                                    '%sstarcat_consensus_clustering.filtered.txt' % self.prefix), '\t')
        spectra_tpm_grouped_hvg.to_csv(os.path.join(self.output_dir, 
                                    '%sstarcat_consensus_spectra_normalized.txt' % self.prefix), '\t')
        spectra_tpm_grouped_hvg.loc[filtsingle_index, :].to_csv(os.path.join(self.output_dir, 
                                    '%sstarcat_consensus_spectra_normalized.filtered.txt' % self.prefix), '\t')
        
        spectra_tpm_grouped.to_csv(os.path.join(self.output_dir, 
                                    '%sstarcat_consensus_spectra_normalized_allgenes.txt' % self.prefix), '\t')
        spectra_tpm_grouped.loc[filtsingle_index, :].to_csv(os.path.join(self.output_dir, 
                                    '%sstarcat_consensus_spectra_normalized_allgenes.filtered.txt' % self.prefix), '\t')        
        
        spectra_scores_grouped.to_csv(os.path.join(self.output_dir, 
                                    '%sstarcat_consensus_spectra_score.txt' % self.prefix), '\t')
        spectra_scores_grouped.loc[filtsingle_index, :].to_csv(os.path.join(self.output_dir, 
                                    '%sstarcat_consensus_spectra_score.filtered.txt' % self.prefix), '\t')
        open(os.path.join(self.output_dir, '%sstarcat_overdispersed_genes_union.txt' % self.prefix), 
             'w').write('\n'.join(hvgs_union))

        top_genes = self.get_top_genes(clus_df, spectra_scores_grouped, n_top_genes=30)
        return(clus_df, spectra_tpm_grouped_hvg, spectra_scores_grouped, hvgs_union, top_genes)
    

    def correlate_geps(self):
        """
        Calculate pairwise correlations of all GEPs using the union highly variable genes between each 2 pairs of cNMF results.

        Parameters
        ----------
        """
        # Calculate all pairwise GEP correlations
        geps = [gep for spectra_tpm in self.spectra_tpm_all for gep in spectra_tpm.index]
        R = pd.DataFrame(np.nan, index = geps, columns = geps)

        for i in range(0, len(self.spectra_tpm_all)):
            for j in range(i, len(self.spectra_tpm_all)):
                # Renormalize each pair of cNMF results using intersecting genes
                overlap_all = list(set(self.spectra_tpm_all[i].columns).intersection(self.spectra_tpm_all[j].columns))

                renorm1 = self.spectra_tpm_all[i][overlap_all].div(
                            self.spectra_tpm_all[i][overlap_all].sum(axis=1), axis=0)*1e6
                renorm2 = self.spectra_tpm_all[j][overlap_all].div(
                            self.spectra_tpm_all[j][overlap_all].sum(axis=1), axis=0)*1e6


                overlap = sorted((set(self.hvgs_all[i]).union(self.hvgs_all[j])).intersection(overlap_all))
        
                renorm1 = renorm1[overlap].div(self.stds_all[i][overlap])
                renorm2 = renorm2[overlap].div(self.stds_all[j][overlap])
                res = df_col_corr(renorm1.T, renorm2.T)
                R.loc[res.index, res.columns] = res
                R.loc[res.columns, res.index] = res.T
        
        self.R = R


    def cluster_geps(self):
        """
        Cluster programs according to the pairwise correlation matrix.

        Parameters
        ----------    

        """

        # Define adjacency matrix between GEP pairs
        self.A = pd.DataFrame(np.zeros((self.R.shape[0], self.R.shape[1])), index = self.R.index, columns = self.R.columns)
        dataset_ind = pd.Series([gep.split(':')[0] for gep in self.R.index], index=self.R.index)
        other_dataset_map = {x:dataset_ind.index[dataset_ind!=x] for x in self.dataset_names}
    
        # Define an edge for correlated GEPs not from the same cNMF result
        for gep in self.R.columns:            
            ds = gep.split(':')[0]
            top_geps = self.R.loc[other_dataset_map[ds], gep].sort_values(ascending = False).head(self.order_thresh)
            top_geps = top_geps.loc[top_geps > self.corr_thresh]

            # Remove GEPs from the same dataset
            tokeep = dataset_ind.loc[top_geps.index].drop_duplicates().index  
            self.A.loc[gep, tokeep] = top_geps.loc[tokeep]

        # Order by correlation
        ord_A = pd.DataFrame(self.A.unstack().sort_values(ascending = False)).reset_index()
        ord_A.columns = ['gep1', 'gep2', 'corr']
        ord_A = ord_A[ord_A['corr']!=0]

        # Cluster GEPs if edge exists between them 
        gep_clusters = pd.DataFrame({'gep' : self.A.index, 'cluster' : 0})
        gep_clusters.index = gep_clusters['gep']

        for i in range(0, ord_A.shape[0]):
            gep = ord_A.loc[i, 'gep1']
            gep_alt = ord_A.loc[i, 'gep2']

            if self.A.loc[gep_alt, gep]>0:
                gep_clusters = self.pairwise_cluster_geps(gep, gep_alt, gep_clusters)

        # Renumber clusters and assign singletons to unique clusters
        clus_dict = {clus_num:sorted(gep_clusters.loc[gep_clusters['cluster']==clus_num, 'gep'].values) 
                 for clus_num in sorted(gep_clusters['cluster'].unique())}
        del clus_dict[0]
        clus_dict_clean = {}
        for clus in clus_dict.keys():
            new_clus = next(i for i, e in enumerate(sorted(clus_dict_clean.keys()) + [ None ], 1) if i != e) 
            clus_dict_clean[new_clus] = clus_dict[clus]
        clus_dict_all = clus_dict_clean.copy()
        for gep in self.R.index:
             if gep not in sum(list(clus_dict_all.values()), []):
                gep_num = next(i for i, e in enumerate(sorted(clus_dict_all.keys()) + [ None ], 1) if i != e) 
                clus_dict_all[gep_num] = [gep]     
        
        # Relabel GEPs and order by cNMF result source
        clus_df = pd.DataFrame.from_dict(clus_dict_all, orient='index', 
                                         columns = ['GEP%d' % x for x in range(1, self.num_results+1)])
        result_names = sorted(clus_df.unstack().dropna().apply(lambda x: x.split(':')[0]).unique())
        
        clus_df_clean = pd.DataFrame(index=clus_df.index, columns=result_names)
        clus_df_clean_forlabel = clus_df.unstack().reset_index().dropna()
        clus_df_clean_forlabel['cNMF_result'] = clus_df_clean_forlabel[0].apply(lambda x: x.split(':')[0])
        for i in clus_df_clean_forlabel.index:
            clus_df_clean.loc[clus_df_clean_forlabel.at[i, 'level_1'], 
                                          clus_df_clean_forlabel.at[i, 'cNMF_result']] = clus_df_clean_forlabel.at[i, 0]
        clus_list = [(clus_df_clean[clus_df_clean==gep]).dropna(how = 'all').index[0] for gep in self.R.index]
        clus_df_clean.index = 'cGEP' + clus_df_clean.index.astype(str)
        clus_df_clean = clus_df_clean[self.dataset_names]

        return(clus_df_clean, clus_list)

    
    def pairwise_cluster_geps(self, gep, gep_alt, gep_clusters):
        '''
        Update clustering given a pair of GEPs. 

        Parameters
        ----------    
        gep: str, name of first program in pair

        gep_alt: str, name of second program in pair

        gep_clusters: pd.DataFrame, dataframe with columns "gep" (program name) and "cluster" (cluster number)

        '''

        # Find which GEPs have been clustered already
        gep_num = gep_clusters.loc[gep, 'cluster']
        alt_gep_num = gep_clusters.loc[gep_alt, 'cluster']
        clustered = [gepX for (gepX, num) in [(gep, gep_num), (gep_alt, alt_gep_num)] if num > 0 ]

        # Merge 2 non-clusterd GEPs
        if len(clustered)==0:
            clus_num = gep_clusters['cluster'].max() + 1
            gep_clusters.loc[[gep, gep_alt], 'cluster'] = clus_num

        # Merge 2 different clusters if enough shared edges and no dataset overlap
        elif len(clustered)==2 and (gep_num != alt_gep_num):
            gep_list = sorted(gep_clusters.loc[gep_clusters['cluster']==gep_num, 'gep'].values) 
            alt_gep_list = sorted(gep_clusters.loc[gep_clusters['cluster']==alt_gep_num, 'gep'].values) 

            num_unconnected = (self.A.loc[gep_list, alt_gep_list]==0).sum().sum() + (
                    self.A.loc[alt_gep_list, gep_list]==0).sum().sum()
            num_total = len(gep_list) * len(alt_gep_list) * 2
            pct_adj = 1 - (num_unconnected / num_total)
            gep_ds = [x.split(':', 1)[0] for x in gep_list]
            alt_gep_ds = [x.split(':', 1)[0] for x in alt_gep_list]

            if ((pct_adj > self.pct_thresh) and (len(set(gep_ds).intersection(alt_gep_ds)) == 0)):
                gep_clusters.loc[alt_gep_list, 'cluster'] = gep_num

        # Add unclustered GEP to a cluster if enough shared edges and no dataset overlap                
        elif len(clustered)==1:
            clus_num = gep_clusters.loc[clustered[0], 'cluster']
            clus_list = sorted(gep_clusters.loc[gep_clusters['cluster']==clus_num, 'gep'].values) 
            unclustered = [gepX for gepX in [gep, gep_alt] if gepX not in clustered]   
            pct_adj = 1 - np.sum((np.array(self.A.loc[unclustered, clus_list] == 0).sum(),
                                  (np.array(self.A.loc[clus_list, unclustered] == 0).sum())))/(len(unclustered)*len(clus_list)*2)
            clus_ds = [x.split(':', 1)[0] for x in clus_list]
            unclustered_ds = [x.split(':', 1)[0] for x in unclustered]
            if (pct_adj > self.pct_thresh) and (len(set(unclustered_ds).intersection(clus_ds)) == 0): 
                gep_clusters.loc[unclustered, 'cluster'] = clus_num   

        elif len(clustered)==2 and (gep_num == alt_gep_num):
            pass
        else:
            sys.exit(-1)

        return(gep_clusters)


    def get_top_genes(self, clus_df, cgep_spectra, n_top_genes=30):
        '''
        Output the top genes for GEPs that got merged into a cGEP and for the cGEP itself

        Parameters
        ----------    
        clus_df: DataFrame, GEP clustering matrix output by cluster_cnmf_results

        cgep_spectra: DataFrame, cGEP spectra scores matrix output by cluster_cnmf_results

        n_top_genes: int, number of top genes to show

        Returns
        ----------  
        top_genes_percgep_dict: dict, {cGEP name : DataFrame of top genes}
        '''
        
        # Get the top genes for each GEP from each dataset separately
        topgenes_perds = {}
        for i,d in enumerate(clus_df.columns):
            topgenes_perds[d] = {}
            top_genes = []
            spectra_scores = self.spectra_score_all[i].T
            
            for gep in spectra_scores.columns:
                top_genes.append(list(spectra_scores.sort_values(by=gep, ascending=False).index[:n_top_genes]))

            topgenes_perds[d] = pd.DataFrame(top_genes, index=spectra_scores.columns, columns=np.arange(1, n_top_genes+1)).T

        # For each cGEP get the top genes from the contributing dataset and the final spectra
        top_genes_percgep_dict = {}
        for cgep in clus_df.index:
            sub_geps = clus_df.loc[cgep, :].dropna()
            top_genes = pd.DataFrame(index=np.arange(1, n_top_genes+1), columns=sub_geps.index)
    
            # Get top genes for each clustered dataset
            for dataset in sub_geps.index:
                top_genes[dataset] = topgenes_perds[dataset][sub_geps.at[dataset]]

            top_genes['cGEP'] = cgep_spectra.loc[cgep, :].sort_values(ascending=False).index[:n_top_genes]
            top_genes_percgep_dict[cgep] = top_genes

        return(top_genes_percgep_dict)