
ANCHOR_COVERAGE_TABLE_DONE?=$(CA_MAP_DIR)/.done_anchor_coverage
$(ANCHOR_COVERAGE_TABLE_DONE): $(CA_ANCHOR_CONTIGS)
	$(_start)
	$(_R) R/ca_info.r anchor.coverage \
		ca.ifn=$(CA_ANCHOR_CONTIGS) \
		cov.ifn=$(COVERAGE_TABLE) \
		ofn=$(ANCHOR_COVERAGE_TABLE)
	$(_end_touch)
anchor_coverage: $(ANCHOR_COVERAGE_TABLE_DONE)

ANCHOR_GC_TABLE_DONE?=$(CA_MAP_DIR)/.done_anchor_gc
$(ANCHOR_GC_TABLE_DONE):
	$(_start)
	$(_R) R/ca_info.r anchor.gc \
		ca.ifn=$(CA_ANCHOR_CONTIGS) \
		binsize=$(ANCHOR_GC_BINSIZE) \
		ifn.gc=$(CONTIG_GC_BINNED) \
		ofn=$(ANCHOR_GC_TABLE)
	$(_end_touch)
anchor_gc: $(ANCHOR_GC_TABLE_DONE)

# median score and median abundance
ANCHOR_PARAMS_DONE?=$(CA_MAP_DIR)/.done_anchor_params
$(ANCHOR_PARAMS_DONE): $(ANCHOR_COVERAGE_TABLE_DONE)
	$(_start)
	$(_R) R/linkage.r anchor.params \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.cov=$(ANCHOR_COVERAGE_TABLE) \
		ofn=$(ANCHOR_PARAMS)
	$(_end_touch)
anchor_params: $(ANCHOR_PARAMS_DONE)

# linkage defined:
# L(c,A) = F(c,A)/median(F(a,A)) * X(A)/X(c)
# where F=O/E and X() is the abundance
LINKAGE_DONE?=$(CA_MAP_DIR)/.done_linkage
$(LINKAGE_DONE): $(ANCHOR_PARAMS_DONE)
	$(_start)
	$(_R) R/linkage.r contig.linkage \
		ifn=$(CA_MATRIX) \
		ifn.anchors=$(ANCHOR_PARAMS) \
		ifn.cov=$(COVERAGE_TABLE) \
		min.anchor.contigs=$(CA_MIN_ANCHOR_CONTIGS) \
		min.contig.coverage=$(CA_CONTIG_COVERAGE) \
		ofn=$(CONTIG_LINKAGE)
	$(_end_touch)
contig_linkage: $(LINKAGE_DONE)

# info:
# - contig median size
# - average gc content
# - number of shared genes
# - relative abundance
# - contact enrichment
# - number of genes
# - number of contigs
ANCHOR_INFO_DONE?=$(CA_MAP_DIR)/.done_info
$(ANCHOR_INFO_DONE):
	$(_start)
	$(_R) R/ca_info.r anchor.info \
		ca.ifn=$(CA_ANCHOR_CONTIGS) \
		ga.ifn=$(CA_ANCHOR_GENES) \
		contig.ifn=$(CONTIG_TABLE) \
		gc.ifn=$(CONTIG_GC_TABLE) \
		coverage.ifn=$(COVERAGE_TABLE) \
		ofn.genes=$(ANCHOR_INFO_TABLE_GENES) \
		ofn.contigs=$(ANCHOR_INFO_TABLE_CONTIGS) \
		ofn.summary=$(ANCHOR_INFO_TABLE)
	$(_end_touch)
ainfo: $(ANCHOR_INFO_DONE)

ANCHOR_INFO_SIZE_DONE?=$(CA_MAP_DIR)/.done_info_size
$(ANCHOR_INFO_SIZE_DONE):
	$(_start)
	$(_R) R/anchor_stats.r anchor.stats.table \
		ifn=$(CA_ANCHOR_CONTIGS) \
		ifn.contigs=$(CONTIG_TABLE) \
		ofn=$(ANCHOR_SIZE_STATS)
	$(_end_touch)
ainfo_size: $(ANCHOR_INFO_SIZE_DONE)

anchor_stats: $(ANCHOR_GC_TABLE_DONE) $(ANCHOR_COVERAGE_TABLE_DONE) $(ANCHOR_INFO_DONE) $(ANCHOR_INFO_SIZE_DONE) $(LINKAGE_DONE)

########################################################################################################################
# plots
########################################################################################################################

plot_anchor_size:
	$(_R) R/plot_basic.r plot.total.size \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.contigs=$(CONTIG_TABLE) \
		fdir=$(ANCHOR_FIGURE_DIR)/total_size

plot_ainfo: $(ANCHOR_INFO_DONE)
	$(_start)
	$(_R) R/ca_info.r plot.anchor.info \
		ifn=$(ANCHOR_INFO_TABLE) \
		ifn.genes=$(ANCHOR_INFO_TABLE_GENES) \
		ifn.contigs=$(ANCHOR_INFO_TABLE_CONTIGS) \
		fdir=$(ANCHOR_FIGURE_DIR)/info
	$(_R) R/ca_info.r plot.anchor.info.summary \
		ifn=$(ANCHOR_INFO_TABLE) \
		fdir=$(ANCHOR_FIGURE_DIR)/info_summary
	$(_end)

plot_linkage:
	$(_start)
	$(_R) R/linkage.r plot.linkage \
		ifn=$(CONTIG_LINKAGE) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		fdir=$(ANCHOR_FIGURE_DIR)/contig_linkage

make_info_plots: plot_anchor_size plot_ainfo plot_linkage

########################################################################################################################
# deprecated
########################################################################################################################

DATASET1?=pre_lib_hic_simple
DATASET2?=post_lib_hic_simple
plot_report:
	$(_R) R/ca_info.r plot.report \
		ifn=$(CA_ANCHOR_CONTIGS) \
		contig.ifn=$(CONTIG_TABLE) \
		gc.ifn=$(CONTIG_GC_TABLE) \
		coverage.ifn1=$(call reval,COVERAGE_TABLE,DATASET=$(DATASET1)) \
		coverage.ifn2=$(call reval,COVERAGE_TABLE,DATASET=$(DATASET2)) \
		gene.ifn=$(GENE_TABLE) \
		uniref.ifn=$(UNIREF_GENE_TAX_TABLE) \
		response.ifn=$(RESPONSE_CONTIG_NORM) \
		contig.cluster.ifn=$(RESPONSE_CONTIG_MEAN_CLUSTER) \
		fdir=$(ANCHOR_FIGURE_DIR)/report

# linkage coeffient of Anchor A is H_A, defined to satisfy:
#  O_c/P_c = 1 + H_A/M_A
ANCHOR_COEFFICENT_DONE?=$(CA_MAP_DIR)/.done_anchor_linkage_coeffient
$(ANCHOR_COEFFICENT_DONE): $(ANCHOR_COVERAGE_TABLE_DONE)
	$(_start)
	$(_R) R/linkage.r anchor.factor \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.cov=$(ANCHOR_COVERAGE_TABLE) \
		ofn=$(ANCHOR_LINKAGE_COEFFICENTS)
	$(_end_touch)
anchor_linkage_coef: $(ANCHOR_COEFFICENT_DONE)

