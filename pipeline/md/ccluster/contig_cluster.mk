#CCLUSTER=$(_md)/bin/ccluster
#CASSIGN=$(_md)/bin/assign_clusters
#CMETRIC=$(_md)/bin/metric

CCLUSTER_BIN?=$(_md)/bin.$(shell hostname)/ccluster
CASSIGN_BIN?=$(_md)/bin.$(shell hostname)/assign_clusters
CMETRIC_BIN?=$(_md)/bin.$(shell hostname)/metric

# obsolete
CMETRIC_DONE?=$(CCLUSTER_DIR)/.done_metric
$(CMETRIC_DONE): $(CONTIG_MATRIX_FILTERED_DONE)
	$(call _start,$(CCLUSTER_DIR))
	$(CMETRIC_BIN) \
	    -contigs $(CONTIG_TABLE) \
            -contacts $(CONTIG_MATRIX_FILTERED) \
            -min_coverage $(CCLUSTER_MIN_COVERAGE) \
            -ratio $(CMETRIC_RATIO) \
            -score_type $(CCLUSTER_SCORE_TYPE) \
            -out $(CCLUSTER_METRIC)
	$(_end_touch)
cmetric: $(CMETRIC_DONE)

CCLUSTER_DONE?=$(CCLUSTER_DIR)/.done_ccluster
$(CCLUSTER_DONE): $(CONTIG_MATRIX_FILTERED_DONE)
	$(call _start,$(CCLUSTER_DIR))
	$(call _time,$(CCLUSTER_DIR)) \
		$(CCLUSTER_BIN) \
		    -contigs $(CONTIG_TABLE) \
	            -contacts $(CONTIG_MATRIX_FILTERED) \
	            -max_contig_pairs $(CCLUSTER_MAX_CONTIG_PAIRS) \
	            -min_score $(CCLUSTER_MIN_SCORE) \
	            -min_coverage $(CCLUSTER_MIN_COVERAGE) \
	            -score_type $(CCLUSTER_SCORE_TYPE) \
		    -cluster_normalize $(CCLUSTER_NORMALIZE_SIZE) \
	            -out_contigs $(CCLUSTER_CONTIGS) \
	            -out_tree $(CCLUSTER_TREE) \
	            -out_scores $(CCLUSTER_SCORES)
	$(_end_touch)
ccluster: $(CCLUSTER_DONE)

##########################################################################################
# assign initial clusters
##########################################################################################

$(INITIAL_ANCHOR_TABLE): $(CCLUSTER_DONE)
	$(call _start,$(CELL_DIR))
	$(CASSIGN_BIN) \
		-contigs $(CCLUSTER_CONTIGS) \
	        -tree $(CCLUSTER_TREE) \
	        -min_score $(CCLUSTER_CUTOFF) \
	        -min_elements $(CCLUSTER_MIN_ELEMENTS) \
	        -out_contigs $@
	$(_end)
make_ianchors: $(INITIAL_ANCHOR_TABLE)
	@echo "DONE contig clustering, INITIAL_ANCHOR_TABLE=$(INITIAL_ANCHOR_TABLE)"

##########################################################################################
# plotting
##########################################################################################

plot_marginals: $(CONTIG_MATRIX_FILTERED)
	$(_start)
	$(_R) R/ccluster.r plot.marginals \
		contigs.ifn=$(CONTIG_TABLE) \
		contacts.ifn=$(CONTIG_MATRIX_FILTERED) \
		threshold=$(CCLUSTER_MIN_COVERAGE) \
		fdir=$(CCLUSTER_FIGURE_DIR)
	$(_end)

plot_gene_graph: $(INITIAL_ANCHOR_TABLE)
	$(_start)
	$(_R) R/ccluster.r plot.complete.graph \
		contigs.ifn=$(CONTIG_TABLE) \
		contacts.ifn=$(CONTIG_MATRIX_FILTERED) \
		fdir=$(CCLUSTER_FIGURE_DIR) \
		anchor.ifn=$(INITIAL_ANCHOR_TABLE) \
		known.ifn=$(CLASSIFY_CONTIG_TABLE)
	$(_end)

plot_cluster_tree: $(INITIAL_ANCHOR_TABLE)
	$(_start)
	$(_R) R/ccluster.r plot.cluster.tree \
		contigs.ifn=$(CCLUSTER_CONTIGS) \
		cluster.table.ifn=$(INITIAL_ANCHOR_TABLE) \
		tree.ifn=$(CCLUSTER_TREE) \
		fdir=$(CCLUSTER_FIGURE_DIR)
	$(_end)


ccluster_init: $(CCLUSTER_BIN) $(CASSIGN_BIN) $(CMETRIC_BIN)

make_plot_metric:
	$(_R) R/contig_matrix.r plot.contig.metric \
		ifn=$(CCLUSTER_METRIC) \
		fdir=$(LIB_FIGURE_DIR)/cc_metric

iplots: plot_cluster_tree make_plot_contig_matrix plot_filter_noise
.PHONY: make_plot_metric plot_cluster_tree plot_gene_graph plot_marginals

$(call _add_module_target,$m,make_ianchors,infer initial anchors)
$(call _add_module_target,$m,iplots,plot initial anchor plots)
