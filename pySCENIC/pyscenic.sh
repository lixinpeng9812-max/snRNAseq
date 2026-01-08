#step1
pyscenic grn --num_workers 12 \
	--sparse \
	--method grnboost2 \
	--output grn.csv \
	sce.loom \
	hs_hgnc_tfs.txt

#step2
pyscenic ctx --num_workers 16 --output regulons.csv --expression_mtx_fname sce.loom --mode "custom_multiprocessing" --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl grn.csv hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

#step3
pyscenic aucell --num_workers 18 --output sample_SCENIC.loom sce.loom regulons.csv
