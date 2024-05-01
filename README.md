# starCAT
Implements *CellAnnotator (aka *CAT/starCAT), annotating scRNA-Seq with predefined gene expression programs

## Citation
If you use *CAT, please cite our [preprint]().

## Installation
To install starCAT, first create a conda environment with necessary dependencies.
```bash
conda create -n cnmf_env --yes --channel bioconda --channel conda-forge --channel defaults python=3.7 fastcluster matplotlib numpy palettable pandas scipy 'scikit-learn>=1.0' pyyaml 'scanpy>=1.8' && conda clean --yes --all # Create environment, cnmf_env, containing required packages
conda activate cnmf_env # Activate cnmf_env - necessary before running cnmf
pip install cnmf # install the actual cnmf package
    
## Only needed to load the example notebook in jupyterlab but not needed for non-interactive runs ## 
conda install --yes jupyterlab && conda clean --yes --all
```

Then install starCAT via the Python Package Index.
```bash
pip install starcatpy
```

## Tutorial
Please see our tutorials in [python](https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vingette.ipynb) and [R](https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vingette_R.ipynb). A sample pipeline using a pre-built reference programs (TCAT.V1) is shown below. 

```python
# Load default TCAT reference from starCAT databse
tcat = starCAT(reference='TCAT.V1')

# tcat.ref.iloc[:5, :5]

#                     A1BG       AARD     AARSD1      ABCA1     ABCB1
# CellCycle-G2M   2.032614  22.965553  17.423538   3.478179  2.297279
# Translation    35.445282   0.000000   9.245893   0.477994  0.000000
# HLA            18.192997  14.632670   2.686475   3.937182  0.000000
# ISG             0.436212   0.000000  18.078197  17.354506  0.000000
# Mito           10.293049   0.000000  52.669895  14.615502  3.341488

# Load cell x genes counts data
adata = tcat.load_counts(datafn)

# Run starCAT
usage, scores = tcat.fit_transform(adata)

usage.iloc[0:2, 0:4]
#                             CellCycle-G2M  Translation       HLA       ISG
# CATGCCTAGTCGATAA-1-gPlexA4       0.000039     0.001042  0.001223  0.000162
# AAGACCTGTAGCGTCC-1-gPlexC6       0.000246     0.100023  0.002991  0.042354

scores.iloc[0:2, :]
#                                  ASA  Proliferation  ASA_binary  \
# CATGCCTAGTCGATAA-1-gPlexA4  0.001556        0.00052       False   
# AAGACCTGTAGCGTCC-1-gPlexC6  0.012503        0.01191       False   

#                             Proliferation_binary Multinomial_Label  
# CATGCCTAGTCGATAA-1-gPlexA4                 False         CD8_TEMRA  
# AAGACCTGTAGCGTCC-1-gPlexC6                 False         CD4_Naive  


```


starCAT also can be run in the command line.
```bash
starcat --reference "TCAT.V1" --counts {counts_fn} --output-dir {output_dir} --name {outuput_name}
```
* --reference - name of a default reference to download (ex. TCAT.V1) OR filepath containing a reference set of GEPs by genes (*.tsv/.csv/.txt), default is 'TCAT.V1'
* --counts - filepath to input (cell x gene) counts matrix as a matrix market (.mtx.gz), tab delimited text file, or anndata file (.h5ad)
* --scores - optional path to yaml file for calculating score add-ons, not necessary for pre-built references
* --output-dir - the output directory. all output will be placed in {output-dir}/{name}...'. default directory is '.'
* --name - the output analysis prefix name, default is 'starCAT'
