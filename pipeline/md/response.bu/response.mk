RES_FDIR?=$(MAP_FIGURE_DIR)/response/$(RESPONSE_ID)

RESPONSE_CONTIG_DONE?=$(RESPONSE_DIR)/.done_contigs
$(RESPONSE_CONTIG_DONE):
	$(call _start,$(RESPONSE_DIR))
	$(_R) R/response.r compute.contig.response \
		ifn=$(CONTIG_TABLE) \
		assembly.dir=$(ASSEMBLY_DIR) \
		ids=$(TDATASETS) \
		ofn.observed=$(RESPONSE_CONTIG_OBSERVED) \
		ofn.expected=$(RESPONSE_CONTIG_EXPECTED) \
		ofn.norm=$(RESPONSE_CONTIG_NORM)
	$(_end_touch)
response_contig: $(RESPONSE_CONTIG_DONE)

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
response: $(RESPONSE_DONE) $(RESPONSE_CONTIG_DONE)

####################################################################
# cutoff contig clustering
####################################################################

$(RESPONSE_CONTIG_TABLE_NORM): $(RESPONSE_CONTIG_DONE)
	$(_R) R/contig_cluster.r norm.conditions \
		ifn=$(RESPONSE_CONTIG_TABLE) \
		ofn=$(RESPONSE_CONTIG_TABLE_NORM)
response_norm: $(RESPONSE_CONTIG_TABLE_NORM)

$(RESPONSE_CONTIG_CLUSTERS): $(RESPONSE_CONTIG_TABLE_NORM)
	$(_start)
	$(call _time,$(RESPONSE_DIR),cluster) $(CORRELATE_CLUSTER_BINARY) \
		-ifn $(RESPONSE_CONTIG_TABLE_NORM) \
		-threshold $(RESPONSE_CLUSTER_THRESHOLD) \
		-method $(RESPONSE_CONTIG_CLUSTER_METHOD) \
		-ofn_bins $(RESPONSE_CLUSTER_HIST) \
		-ofn $@

cluster_contigs: $(RESPONSE_CONTIG_CLUSTERS)

FILES_R=$(addprefix $(_md)/cpp/,Params.cpp Params.h util.h util.cpp)
$(eval $(call bin_rule2,correlation_cluster,$(FILES_R)))

CORRELATE_CLUSTER_BINARY=$(_md)/bin/correlation_cluster

####################################################################
# kmeans contig clustering
####################################################################

$(eval $(call bin_rule2,coverage_kmeans,$(FILES_R)))
COVERAGE_KMEANS_BINARY=$(_md)/bin/coverage_kmeans

# convert to fold change
$(RESPONSE_CONTIG_TABLE_FOLD): $(RESPONSE_CONTIG_DONE)
	$(_R) R/contig_cluster.r fold.change \
		ifn=$(RESPONSE_CONTIG_TABLE) \
		max.change=$(RESPONSE_KMEANS_MAX_FOLD) \
		min.coverage=$(RESPONSE_KMEANS_MIN_CONTIG_COVERAGE) \
		ofn=$(RESPONSE_CONTIG_TABLE_FOLD)
response_fold: $(RESPONSE_CONTIG_TABLE_FOLD)

# cluster response using kmeans
RESPONSE_KMEANS_DONE?=$(RESPONSE_DIR)/.done_kmeans
$(RESPONSE_KMEANS_DONE): $(RESPONSE_CONTIG_TABLE_FOLD)
	$(_start)
	$(call _time,$(RESPONSE_DIR),cluster) $(COVERAGE_KMEANS_BINARY) \
		-ifn $(RESPONSE_CONTIG_TABLE_FOLD) \
		-n_clusters $(RESPONSE_KMEANS_N_CLUSTERS) \
		-random_seed $(RESPONSE_KMEANS_SEED) \
		-cluster_seed_method $(RESPONSE_KMEANS_SEED_METHOD) \
		-ofn_prefix $(RESPONSE_KMEANS_PREFIX)
	$(_end_touch)
response_kmeans: $(RESPONSE_KMEANS_DONE)

$(RESPONSE_KMEANS_ORDER): $(RESPONSE_KMEANS_DONE)
	$(_R) R/contig_cluster.r cluster.order \
		ifn=$(RESPONSE_KMEANS_CENTROID_MEAN) \
		ofn=$@
response_cluster_order: $(RESPONSE_KMEANS_ORDER)

####################################################################
# plots
####################################################################

plot_response_kmeans:
	$(_R) R/contig_cluster.r plot.centroids.colors \
		ifn.mean=$(RESPONSE_KMEANS_CENTROID_MEAN) \
		ifn.contigs=$(RESPONSE_KMEANS_CONTIGS) \
		ifn.contig.table=$(CONTIG_TABLE) \
		ifn.order=$(RESPONSE_KMEANS_ORDER) \
		fdir=$(RES_FDIR)/kmeans


plot_response: $(RESPONSE_DONE)
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

