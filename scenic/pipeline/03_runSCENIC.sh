#!/bin/bash

cd ~/test_scenic_metabolic

echo "pyscenic ctx"

pyscenic ctx Output/W.regular-metabolite.adj.tsv \
  ../2023-12-08_metabolicSCENIC/auxilliaries/dm6-5kb-upstream-full-tx-11species.mc8nr.genes_vs_motifs.rankings.feather \
  --annotations_fname ../2023-12-08_metabolicSCENIC/resources/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl \
  --expression_mtx_fname Data/8W_exprMat_regular.loom \
  --output Output/W.regulons.tsv --num_workers 12

pyscenic ctx Output/Y.regular-metabolite.adj.tsv \
  ../2023-12-08_metabolicSCENIC/auxilliaries/dm6-5kb-upstream-full-tx-11species.mc8nr.genes_vs_motifs.rankings.feather \
  --annotations_fname ../2023-12-08_metabolicSCENIC/resources/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl \
  --expression_mtx_fname Data/8Y_exprMat_regular.loom \
  --output Output/Y.regulons.tsv --num_workers 12

echo "pyscenic aucell"

pyscenic aucell \
  Data/8W_exprMat_regular.loom \
  Output/W.regulons.tsv \
  -o Output/8W.auc_mtx.loom --num_workers 12

pyscenic aucell \
  Data/8Y_exprMat_regular.loom \
  Output/Y.regulons.tsv \
  -o Output/8Y.auc_mtx.loom --num_workers 12
