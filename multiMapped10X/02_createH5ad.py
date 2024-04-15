# %%
import pandas as pd
from scipy.sparse import csr_matrix
import anndata
import os
import gzip
import re

barcode_file="barcodes.tsv"
feature_file="features.tsv"

barcodes = pd.read_csv(barcode_file, sep='\t', header=None).iloc[:,0]
features = pd.read_csv(feature_file, sep='\t', header=None).iloc[:,0]

def load_sparse_matrix(file_path):
    rows = []
    cols = []
    data = []

    with gzip.open(file_path, 'rt') as file:
        
        for line in file:
            if line.startswith('%'):
                continue
            
            parts = line.strip().split()
            row = int(parts[0])
            col = int(parts[1])
            value = float(parts[2])
            rows.append(row)
            cols.append(col)
            data.append(value)
    
    # get the first row(summary) out 
    total_read = data.pop(0)
    n_rows = rows.pop(0)
    n_cols = cols.pop(0)

    # print h5 file size
    print("The first row of mtx file are : " + str(n_rows) + "," + str(n_cols)+ "," + str(total_read))       

    # Create csr format sparse matrix  
    # add 1 as python is 0-based indexed   
    sparse_matrix = csr_matrix((data, (rows, cols)), shape=(n_rows+1, n_cols+1))

    return sparse_matrix

file_path = "UniqueAndMult-EM.mtx.gz"
sparse_matrix = load_sparse_matrix(file_path)

# Create AnnData object
adata = anndata.AnnData(sparse_matrix)
print(adata.shape)

# Write to h5ad file
adata.write('UniqueAndMult.h5ad')
