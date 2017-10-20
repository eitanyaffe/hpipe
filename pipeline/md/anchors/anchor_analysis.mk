pshared:
	$(_R) R/shared_genes.r plot.shared.genes \
		ca.matrix.ifn=$(CA_MATRIX) \
		contigs.ifn=$(CA_ANCHOR_CONTIGS) \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		genes.ifn=$(GENE_TABLE) \
		uniref.ifn=$(UNIREF_GENE_TAX_TABLE)
