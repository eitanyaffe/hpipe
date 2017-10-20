# anchor snp density
plot_anchor_snp:
	$(_R) $(_md)/R/anchor_variance.r plot.snp.density \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		assembly.dir=$(ASSEMBLY_DIR) \
		dataset1=$(VARIANCE_DATASET1) \
		dataset2=$(VARIANCE_DATASET2) \
		fdir=$(CA_MAP_FDIR)/anchor_variance
