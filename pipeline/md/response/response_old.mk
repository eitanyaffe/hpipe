####################################################################
# contig coverage matrix
####################################################################

RESPONSE_CONTIG_DONE?=$(RESPONSE_DIR)/.done_contigs
$(RESPONSE_CONTIG_DONE):
	$(call _start,$(RESPONSE_DIR))
	$(_R) R/response.r compute.contig.response \
		ifn=$(CONTIG_TABLE) \
		do.log=$(RESPONSE_LOG) \
		reg=$(RESPONSE_REG) \
		assembly.dir=$(ASSEMBLY_DIR) \
		ids=$(TDATASETS) \
		ofn.observed=$(RESPONSE_CONTIG_OBSERVED) \
		ofn.expected=$(RESPONSE_CONTIG_EXPECTED) \
		ofn.norm=$(RESPONSE_CONTIG_NORM)
	$(_end_touch)
response_contig: $(RESPONSE_CONTIG_DONE)

####################################################################
# contig clustering
####################################################################

RESPONSE_MEAN_CLUSTER_DONE?=$(RESPONSE_DIR)/.done_mean_cluster
$(RESPONSE_MEAN_CLUSTER_DONE): $(RESPONSE_CONTIG_DONE)
	$(call _start,$(RESPONSE_DIR))
	$(_R) R/response.r mean.cluster \
		ifn=$(RESPONSE_CONTIG_NORM) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		thresholds=$(RESPONSE_CLUSTER_THRESHOLDS) \
		ofn.prefix=$(RESPONSE_CONTIG_MEAN_CLUSTER_PREFIX) \
		ofn.order=$(RESPONSE_CONTIG_ORDER)
	$(_end_touch)
cluster_contigs: $(RESPONSE_MEAN_CLUSTER_DONE)

RESPONSE_INFO_DONE?=$(RESPONSE_SDIR)/.done_info
$(RESPONSE_INFO_DONE): $(RESPONSE_MEAN_CLUSTER_DONE)
	$(call _start,$(RESPONSE_SDIR))
	$(_R) R/response_info.r info \
		clusters.ifn=$(RESPONSE_CONTIG_MEAN_CLUSTER) \
		contig.ifn=$(CONTIG_TABLE) \
		coverage.ifn=$(COVERAGE_TABLE) \
		gc.ifn=$(CONTIG_GC_TABLE) \
		gene.ifn=$(GENE_TABLE) \
		uniref.ifn=$(UNIREF_GENE_TAX_TABLE) \
		ca.ifn=$(CA_ANCHOR_CONTIGS) \
		ofn=$(RESPONSE_INFO)
	$(_end_touch)
response_info: $(RESPONSE_INFO_DONE)

RESPONSE_PATTERN_DONE?=$(RESPONSE_SDIR)/.done_pattern
$(RESPONSE_PATTERN_DONE): $(RESPONSE_MEAN_CLUSTER_DONE)
	$(call _start,$(RESPONSE_SDIR))
	$(_R) R/response.r pattern.response \
		ifn.contigs=$(RESPONSE_CONTIG_NORM) \
		ifn.clusters=$(RESPONSE_CONTIG_MEAN_CLUSTER) \
		ofn.mean=$(RESPONSE_PATTERN_MEAN) \
		ofn.sd=$(RESPONSE_PATTERN_SD)
	$(_end_touch)
pattern: $(RESPONSE_PATTERN_DONE)

# match clusters and anchors
RESPONSE_ANCHOR_DONE?=$(RESPONSE_SDIR)/.done_anchor
$(RESPONSE_ANCHOR_DONE): $(RESPONSE_PATTERN_DONE)
	$(_start)
	$(_R) R/response.r anchor.pattern \
		ifn.contigs=$(CONTIG_TABLE) \
		ifn.clusters=$(RESPONSE_CONTIG_MEAN_CLUSTER) \
		ifn.anchors=$(CA_ANCHOR_CONTIGS) \
		ofn=$(RESPONSE_EA_MATRIX)
	$(_end_touch)
anchor_pattern: $(RESPONSE_ANCHOR_DONE)

response: $(RESPONSE_INFO_DONE) $(RESPONSE_ANCHOR_DONE)

####################################################################
# anchor response (TBD)
####################################################################

RESPONSE_DONE?=$(RESPONSE_DIR)/.done
$(RESPONSE_DONE):
	$(call _start,$(RESPONSE_DIR))
	$(_R) R/temporal_plots.r compute.response \
		cc.table=$(ANCHOR_TABLE) \
		assembly.dir=$(ASSEMBLY_DIR) \
		ids=$(TDATASETS) \
		ofn.response=$(RESPONSE_TABLE) \
		ofn.order=$(RESPONSE_ORDER)
	$(_end_touch)
	@mkdir -p $(ORDER_DIR)
	cp -r $(RESPONSE_ORDER) $(ORDER_DIR)/response_$(RESPONSE_ID)
response_old: $(RESPONSE_DONE) $(RESPONSE_CONTIG_DONE)

####################################################################
# plots
####################################################################

RES_FDIR?=$(MAP_FIGURE_DIR)/response/$(RESPONSE_ID)

plot_patterns:
	$(_R) R/response_plot.r plot.pattern.lines \
		ifn.mean=$(RESPONSE_PATTERN_MEAN) \
		ifn.sd=$(RESPONSE_PATTERN_SD) \
		ifn.clusters=$(RESPONSE_CONTIG_CLUSTERS) \
		fdir=$(RES_FDIR)/patterns

plot_response_kmeans:
	$(_R) R/contig_cluster.r plot.centroids.colors \
		ifn.mean=$(RESPONSE_KMEANS_CENTROID_MEAN) \
		ifn.contigs=$(RESPONSE_KMEANS_CONTIGS) \
		ifn.contig.table=$(CONTIG_TABLE) \
		ifn.order=$(RESPONSE_KMEANS_ORDER) \
		fdir=$(RES_FDIR)/kmeans


plot_anchor_response: $(RESPONSE_DONE)
	$(_R) R/temporal_plots.r plot.response \
		ids=$(TDATASETS) \
		ifn.response=$(RESPONSE_TABLE) \
		taxa.ifn=$(SET_TAXA_REP_LEGEND) \
		order.ifn=$(RESPONSE_ORDER) \
		fdir=$(RES_FDIR)

plot_time_class: $(CC_GENE_TAX_TABLE)
	$(_R) R/temporal_plots.r plot.temporal.genes \
		cc.table=$(CC_MAP_TABLE) \
		gene.table=$(CC_GENE_TAX_TABLE) \
		tdatasets=$(TDATASETS) \
		assembly.dir=$(ASSEMBLY_DIR) \
		odir=$(RES_FDIR)

response_plots: plot_response

response_init: $(CORRELATE_CLUSTER_BINARY) $(COVERAGE_KMEANS_BINARY)

