# starCAT
Implements *CellAnnotator (aka *CAT), annotating scRNA-Seq with predefined gene expression programs.

## Citation
If you use *CAT, please cite our preprint.

## Tutorial
Please see our tutorials in [python](https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vingette.ipynb) and [R](https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vingette_R.ipynb). Briefly, a sample pipeline using a pre-built reference programs (TCAT.V1) is shown below. 

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
