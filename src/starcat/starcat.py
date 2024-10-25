#!/usr/bin/env python

import pandas as pd
import numpy as np
from anndata import read_h5ad, AnnData
import os
import scipy.sparse as sp
import yaml
import requests
import warnings
from sklearn.decomposition import non_negative_factorization
from sklearn import preprocessing
from .utils import read_10x_mtx, decompress_tar

reference_url = os.path.join(os.path.dirname(__file__), 'current_references.tsv')
ref_list = pd.read_csv(reference_url, comment='#', sep='\t')
available_refs = ref_list['Name'].values

_nmf_kwargs = dict(
                   beta_loss='frobenius',
                   solver='mu',
                   tol=1e-4,
                   max_iter=1000,
                   init='random',
                   update_H = False
                   )

def is_integer_matrix(matrix):
    # Check if the matrix is a sparse matrix
    if sp.issparse(matrix):
        # Convert to a dense format without changing the data
        data = matrix.data
    else:
        # It's a numpy array
        data = np.ravel(matrix)  # Flatten the matrix to a 1D array for easy processing
    
    # Check each element if it's an integer or a float equivalent to an integer
    return np.all(data == np.floor(data))


def is_nonnegative_matrix(matrix):
    # Check if the matrix is a sparse matrix
    if sp.issparse(matrix):
        # Convert to a dense format without changing the data
        data = matrix.data
    else:
        # It's a numpy array
        data = np.ravel(matrix)  # Flatten the matrix to a 1D array for easy processing
    
    # Check each element if it's an integer or a float equivalent to an integer
    return np.all(data >= 0)


def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj


class starCAT():  
    def __init__(self, reference = 'TCAT.V1', score_path = None, cachedir='./cache'):
        """
        Runs *CAT on a query dataset.
        
        Parameters
        ----------
        reference : str, name of reference from starCAT database or path to a custom reference file ending in .tsv or .txt (Default="TCAT.V1") 
        
        score_path : str, optional path to a yaml file to compute addon scores for discrete or continuous features from starCAT output.
        Only needed if reference is a path to a custom reference file and not a starCAT database name (default=None)        
        """
        
        self._nmf_kwargs = _nmf_kwargs
        self._cache = cachedir
        self.ref_name = reference
        self.usage = None
        self.usage_norm = None
        self.scores = None
        
        if reference.endswith('.txt') or reference.endswith('.tsv'):
            print("Using user specified reference spectra file {f}".format(f=reference))
            self.ref = pd.read_csv(reference, index_col = 0, sep = '\t').astype(np.float32)
            self.score_data = self.load_scores(score_path)
            self.score_path = score_path
        else:
            # Download cache if necessary and load references and scores from cache
            print("Using reference from starCAT database")
            self.ref, self.score_data, self.score_path = self._initialize_ref()
        
        #Filter genes with NaN in them from reference
        num_nulls = self.ref.isnull().sum(axis=0)
        has_null = num_nulls>0
        nnull_gs = has_null.sum()
        nref_gs = self.ref.shape[1]
        if nnull_gs > 0:
            warnings.warn("Filtering %d of %d reference genes that contain NaN for 1+ GEPs." % (nnull_gs, nref_gs))
            self.ref = self.ref.loc[:,~has_null]


    def _initialize_ref(self):
        """"
        Initializes a reference set of GEPs by genes from a file name ending in the extension .txt or .tsv
        or using the name of a pre-built reference (default reference is TCAT.V1).
        
        """

        if self.ref_name not in available_refs:
            refliststr = ','.join(available_refs)
            raise Exception(
                "{ref} is not found in list of pre-built reference names. "
                "It is also not a valid path to a reference file which would "
                "need to end in .tsv or .txt. Please provide a valid file path "
                "or a reference string from among the following [{posrefs}]".format(ref=self.ref_name,
                                                                                    posrefs=refliststr)
            )

        spectra_fn = os.path.join(self._cache, self.ref_name, self.ref_name + '.reference.tsv')
        score_path = os.path.join(self._cache, self.ref_name, '%s.scores.yaml' % self.ref_name)

        if os.path.exists(spectra_fn):
            print('Loading reference from existing cache file for reference %s' % self.ref_name)
        else:
            print('Downloading reference %s to cache' % self.ref_name)
            ref_url = ref_list.loc[ref_list['Name']==self.ref_name, 'Link'].values[0]
            self._download_to_cache(ref_url)
       
        ref = pd.read_csv(spectra_fn, index_col = 0, sep = '\t').astype(np.float32)

        if os.path.exists(score_path):
            scoredat = self.load_scores(score_path)
        else:
            print('No derived score data available for reference %s' % self.ref_name)
            scoredat = None
            
        return(ref, scoredat, score_path)


    def load_scores(self, score_path = None):
        '''
        Loads score functions from a .yaml file.
        '''
        if score_path == None:
            print('No scores provided')
            return(None)
        else:
            with open(score_path, 'r') as f:
                return(yaml.safe_load(f))
        

    def _download_to_cache(self, ref_url):
        # Cache reference 
        if not os.path.exists(self._cache):
            os.makedirs(self._cache)
            print('Making empty cache directory "%s"' % self._cache)

        if not os.path.isdir(os.path.join(self._cache, self.ref_name)):
            response = requests.get(ref_url, stream=True)
            response.raise_for_status()
                
            print('Caching reference to %s' % os.path.join(self._cache, self.ref_name))
            tar_fn = os.path.join(self._cache, '%s.tar.gz' % self.ref_name)
            
            with open(tar_fn, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)    

            decompress_tar(tar_fn, outdir=self._cache)
            os.remove(tar_fn)
    
    def load_counts(self, counts_fn):
        """
        Loads counts matrix (cells X genes) from 10X outputs, tab delimited text file, or anndata file.
        
        Parameters:
        ----------
        counts_fn : str, path to input counts matrix (ex. *.h5ad, *.mtx.gz). If neither, assumes tab-delimited text file
        
        """
        if not os.path.exists(counts_fn):
            raise FileNotFoundError(f"The specified file does not exist: {counts_fn}")

        if counts_fn.endswith('.h5ad'):
            adata = read_h5ad(counts_fn) 
            
        elif counts_fn.endswith('.mtx.gz') or counts_fn.endswith('.mtx'):
            counts_dir = os.path.dirname(counts_fn)
            adata = read_10x_mtx(path = counts_dir)
                
        # Convert other forms of query data to AnnData objects
        else:
            input_counts = pd.read_csv(counts_fn, sep='\t', index_col=0)
            adata = AnnData(X=input_counts.values,
                                   obs=pd.DataFrame(index=input_counts.index),
                                   var=pd.DataFrame(index=input_counts.columns))
            
        return adata

                        
    def fit_transform(self, query, return_unnormalized=False):
        """
        Takes an input data matrix and a fixed spectra and uses NNLS to find the optimal
        usage matrix. If input data are pandas.DataFrame, returns a DataFrame with row
        index matching X and columns index matching index of spectra

        Parameters
        ----------
        query : query AnnData object or pandas.DataFrame, cells X genes
        
        return_unnormalized : boolean, if True, return unnormalized usage. Default returns usages normalized to 1

        """

        query = self.prepare_query(query)        
        self.usage = self.fit_query_usage(query)
        self.usage_norm = self.usage.div(self.usage.sum(axis=1), axis=0)
        
        if self.score_data is not None:
            self.scores = self.compute_addon_scores()
        
        if return_unnormalized:
            return self.usage, self.scores
        else:
            return self.usage_norm, self.scores
            
    
    
    def prepare_query(self, query):
        """
        Load query dataset as AnnData object and optionally perform normalization.

        adata : query AnnData object

        """
        
        if not isinstance(query, (pd.DataFrame, AnnData)):
            raise TypeError('%s is not a valid object type. Please convert to pd.DataFrame or AnnData.' % type(query))
        elif isinstance(query, pd.DataFrame):
            query = AnnData(X=query.values, obs=pd.DataFrame(index=query.index), var=pd.DataFrame(index=query.columns))
            
        if not is_nonnegative_matrix(query.X):
            raise Exception("""query input contains negative values. Must be non-negative""")
            
        if not is_integer_matrix(query.X):
            warnings.warn("""WARNING!: query input is not an integer count matrix as expected.
            Please provide an integer count matrix unless you are sure you know what you are doing.""", UserWarning)
            
        overlap_genes = sorted(set(self.ref.columns).intersection(set(query.var.index)))
        print('%d out of %d genes in the reference overlap with the query' % (len(overlap_genes), self.ref.shape[1]))
        self.overlap_genes = overlap_genes
        query = query[:, self.overlap_genes].copy()
        num_zeros = (np.array(query.X.sum(axis=0)).reshape(-1) == 0).sum()
        if num_zeros > 0:
            warnings.warn("""WARNING!: query input has %d genes with 0 counts after overlapping with query. Normalized values for these genes are set to 0.""" % num_zeros, UserWarning)
        
        query.X = preprocessing.scale(query.X, with_mean=False)

        # In python 3.10, scale casts query as float64. Here just ensure its same dtype as reference
        if query.X.dtype != self.ref.values.dtype:
            query.X = query.X.astype(self.ref.values.dtype)
        
        return query
         
        
    def fit_query_usage(self, query):
        rf_usages = self.refit_usage(query.X, self.ref[self.overlap_genes].values,
                         self._nmf_kwargs.copy())          
        rf_usages = pd.DataFrame(rf_usages, index=query.obs.index,
                                 columns=self.ref.index)
        return(rf_usages)
        
        
    def refit_usage(self, X, spectra, nmf_kwargs):
        """
        Takes an input data matrix and a fixed spectra and uses NNLS to find the optimal
        usage matrix. If input data are pandas.DataFrame, returns a DataFrame with row
        index matching X and columns index matching index of spectra

        Parameters
        ----------
        X : pandas.DataFrame or numpy.ndarray, cells X genes
            Non-negative expression data to fit spectra to

        spectra : pandas.DataFrame or numpy.ndarray, programs X genes
            Non-negative spectra of expression programs
        """

        nmf_kwargs['H'] = spectra
        nmf_kwargs['n_components'] = spectra.shape[0]
        rf_usages, _, _ = non_negative_factorization(X, **nmf_kwargs)
        return(rf_usages)
    
    
    def compute_addon_scores(self):
        """
        Compute pre-built reference score add-ons to usages. Currently can be applied to pre-built reference datasets only.
        """

        if self.score_data is None:
            raise Exception("""No addon scores were provided. Please provide a score_path input or use a default
            reference that includes addon scores when initializing starCAT""")
            
        if self.usage is None:
            raise Exception("""No usages are available. Please run `fit_transform` first to fit usages.""")            
                
        score_res = pd.DataFrame()
        for score_type in self.score_data['scores']:
            for score in self.score_data['scores'][score_type]:                
                if score['normalization']=='normalized':
                    usage_for_score = self.usage_norm
                else:
                    usage_for_score = self.usage

                # By default define each score as a weighted sum of GEPs, with default weights of 1
                if 'file' not in score.keys():
                    if 'weights' in score.keys():
                        weights = score['weights']
                    else:
                        weights = np.array([1]*len(score['columns']))
                    score_values = pd.Series(np.dot(usage_for_score[score['columns']], weights), 
                                             index = usage_for_score.index, name = score['name'])

                    if score_type == 'continuous':
                        score_res[score['name']] = score_values
                    else:
                        score_res[score['name']] = score_values > score['threshold']    
                        
                # Optionally execute score code blocks
                else:
                    global_vars = {}
                    local_vars = {}
                    global_vars['usage_for_score'] = usage_for_score
                    global_vars['score_dir'] = os.path.dirname(self.score_path) 
                    
                    with open(os.path.join(global_vars['score_dir'], score['file']), 'r') as F:
                        exec(F.read(), global_vars, local_vars)
                    score_res[score['name']] = local_vars['res']
        
        return score_res
    
    
    def save_results(self, output_dir, name):
        """
        Output *CAT usages and add-on scores to a txt file located in [output-dir]/[name]...
        """
        # self.usage.to_csv(os.path.join(output_dir, name+'.rf_usage.txt'), sep='\t')
        self.usage_norm.to_csv(os.path.join(output_dir, name+'.rf_usage_normalized.txt'), sep='\t')
        print('Saving usages to %s' % os.path.join(output_dir, name+'.rf_usage_normalized.txt'))
        
        if not self.scores is None:
            self.scores.to_csv(os.path.join(output_dir, name+'.scores.txt'), sep='\t')
            print('Saving scores to %s' % os.path.join(output_dir, name+'.scores.txt'))


def main():
    """
    Example command:
    
    python starCAT.py --reference "TCAT.V1" --counts $counts --output-dir $output_dir --name $name
    
    """
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', type=str, help='File containing a reference set of GEPs by genes (*.tsv/.csv/.txt) OR the name of a default reference to use (ex. TCAT.V1).', default='TCAT.V1')
    parser.add_argument('-c', '--counts', type=str, help='Input (cell x gene) counts matrix as df.npz, tab delimited text file, or anndata file (h5ad)')
    parser.add_argument('--output-dir', type=str, help='Output directory. All output will be placed in [output-dir]/[name]...', nargs='?', default='.')
    parser.add_argument('--name', type=str, help='Name for analysis. All output will be placed in [output-dir]/[name]...', nargs='?', default='starCAT')    
    parser.add_argument('-s', '--scores', type=str, help='Optional path to yaml file for calculating score add-ons. Not necessary for pre-built references', default=None)

        
    args = parser.parse_args()
    cat_obj = starCAT(reference = args.reference, score_path = args.scores)
    adata = cat_obj.load_counts(counts_fn = args.counts)
    rf_usage, score_res = cat_obj.fit_transform(query = adata)
    cat_obj.save_results(args.output_dir, args.name)
    
    
if __name__=="__main__":
    main()
