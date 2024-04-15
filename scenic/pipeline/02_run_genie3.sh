#!/bin/bash

cd ~/test_scenic_metabolic
mkdir -p Output
echo "pyscenic genie3"

# run genie3 as separate step 
pyscenic grn --num_workers 12 -m genie3 \
  -o Output/Y.regular-metabolite.adj.tsv \
  Data/8Y_exprMat_regular.loom \
  resources/TFs/metabolic_TFs_dmel.txt

pyscenic grn --num_workers 12 -m genie3 \
  -o Output/W.regular-metabolite.adj.tsv \
  Data/8W_exprMat_regular.loom \
  resources/TFs/metabolic_TFs_dmel.txt

