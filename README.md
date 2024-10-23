## starCAT <img src="https://drive.google.com/uc?export=view&id=1W1in9vldkKdNe6ncwsHD6L6MSvfcKV6M" width="130px" align="right" />
Implements *CellAnnotator (aka *CAT/starCAT), annotating scRNA-Seq with predefined gene expression programs
<br>

## Citation

If you use *CAT, please cite our [preprint](https://doi.org/10.1101/2024.05.03.592310).

## Installation

You can install starCAT and its dependencies via the Python Package Index.
```bash
pip install starcatpy
```

We tested it with scikit-learn 1.3.2, AnnData 0.9.2, and python 3.8. To run the tutorials, you also need jupyter or jupyterlab as well as scanpy and cnmf:

```bash
pip install jupyterlab scanpy cnmf
```


## Basic starCAT usage
Please see our tutorials in [python](Examples/starCAT_vignette.ipynb) and [R](Examples/starCAT_vignette_R.ipynb). A sample pipeline using a pre-built reference programs (TCAT.V1) is shown below. 

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

## Alternate implementation
For small datasets (<2 GB), try running starCAT without installing any packages on our [website](https://immunogenomics.io/starcat/).

## Creating your own reference

We provide example scripts for constructing custom starCAT references from [a single cNMF run](./Examples/build_reference_vignette.ipynb) or [multiple cNMF runs](./Examples/build_multidataset_reference_vignette.ipynb). 

__Please let us know if you are interested in making your reference publically available for others to use analogous to our TCAT.V1 reference. You can email me at dkotliar@broadinstitute.org__
