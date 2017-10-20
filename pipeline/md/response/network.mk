################################################################################
# element-anchor association
################################################################################

#EE_MAP_DONE?=$(EE_BASE_MAP_DIR)/.done_ee_map
#$(EE_MAP_DONE):  $(RESPONSE_CLUSTER_DONE)
#	$(_start)
#	$(MAKE) m=anchors ee_map LIB_ID=$(RESPONSE_LIB_ID) EE_QSUB_DIR=$(QSUB_DATASET_DIR)/temporal/$(RESPONSE_LIB_ID)
#	$(_end_touch)
#response_map: $(EE_MAP_DONE)

response_matrix_internal:
	$(call _start,$(EE_MAP_DIR))
	$(_R) R/merge_ca.r merge.ca \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.contigs=$(CONTIG_TABLE) \
		ifn.clusters=$(RESPONSE_CONTIG_CLUSTER) \
		ifn.ca=$(CA_MATRIX) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		min.enrichment=$(RESPONSE_NETWORK_MIN_ENRICHMENT) \
		separate.min.contacts=$(RESPONSE_NETWORK_SEPARATE_MIN_CONTACTS) \
		separate.max.enrichment=$(RESPONSE_NETWORK_SEPARATE_MAX_SCORE) \
		ofn=$(RESPONSE_ANCHOR_MATRIX)
	$(_end)

RESPONSE_MATRIX_DONE?=$(EE_BASE_MAP_DIR)/.done_ca_merge
$(RESPONSE_MATRIX_DONE):
	$(_start)
	$(MAKE) response_matrix_internal LIB_ID=$(RESPONSE_LIB_ID) ANCHOR_CLUSTER_TABLE=$(ANCHOR_CLUSTER_TABLE)
	$(_end_touch)
response_matrix: $(RESPONSE_MATRIX_DONE)

RESPONSE_ANCHOR_ELEMENTS_DONE?=$(EE_BASE_MAP_DIR)/.done_anchor_elements
$(RESPONSE_ANCHOR_ELEMENTS_DONE): $(RESPONSE_MATRIX_DONE)
	$(_start)
	$(_R) R/response.r anchor.elements \
		ifn.matrix=$(RESPONSE_ANCHOR_MATRIX) \
		ifn.means=$(RESPONSE_PATTERN_MEAN) \
		ofn=$(RESPONSE_ANCHOR_ELEMENTS)
	$(_end_touch)
response_elements: $(RESPONSE_ANCHOR_ELEMENTS_DONE)

# all network
RESPONSE_NETWORK_DONE?=$(EE_BASE_MAP_DIR)/.done_network
$(RESPONSE_NETWORK_DONE): $(RESPONSE_ANCHOR_ELEMENTS_DONE)
	$(_start)
	$(_R) R/network.r compute.network \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_ANCHOR_ELEMENTS) \
		ifn.contigs=$(RESPONSE_CONTIG_CLUSTER) \
		ifn.map=$(RESPONSE_ANCHOR_MATRIX) \
		ifn.element.info=$(RESPONSE_INFO) \
		ifn.identity=$(ANCHOR_MATRIX_TABLE) \
		ifn.genes=$(RESPONSE_GENES) \
		ofn.network=$(RESPONSE_NETWORK) \
		ofn.elements=$(RESPONSE_NETWORK_ELEMENTS) \
		ofn.genes=$(RESPONSE_NETWORK_GENES)
	$(_end_touch)
network: $(RESPONSE_NETWORK_DONE)

# limit to (a) >1 host (b) minimal number of uniref genes
RESPONSE_SELECT_NETWORK_DONE?=$(EE_BASE_MAP_DIR)/.done_network_select
$(RESPONSE_SELECT_NETWORK_DONE): $(RESPONSE_NETWORK_DONE)
	$(_start)
	$(_R) R/network.r select.network \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.network=$(RESPONSE_NETWORK) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS) \
		ifn.genes=$(RESPONSE_NETWORK_GENES) \
		min.hosts=$(NETWORK_MIN_HOSTS) \
		min.genes=$(NETWORK_MIN_GENES) \
		min.uniref.genes=$(NETWORK_MIN_UNIREF_GENES) \
		max.pearson=$(NETWORK_MAX_PEARSON) \
		max.diameter=$(NETWORK_MAX_DIAMETER) \
		ofn.network=$(RESPONSE_NETWORK_SELECT) \
		ofn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ofn.genes=$(RESPONSE_NETWORK_GENES_SELECT)
	$(_end_touch)
network_select: $(RESPONSE_SELECT_NETWORK_DONE)

RESPONSE_COORD_NETWORK_DONE?=$(EE_BASE_MAP_DIR)/.done_network_coord
$(RESPONSE_COORD_NETWORK_DONE): $(RESPONSE_SELECT_NETWORK_DONE)
	$(_start)
	$(_R) R/network.r network.coord \
		ifn=$(RESPONSE_NETWORK_SELECT) \
		type=$(RESPONSE_NETWORK_LAYOUT_TYPE) \
		anchor.order.ifn=$(ANCHOR_CLUSTER_TABLE) \
		ofn=$(RESPONSE_NETWORK_COORDS)
	$(_end_touch)
network_coord: $(RESPONSE_COORD_NETWORK_DONE)

RESPONSE_GENE_STATS_DONE?=$(EE_BASE_MAP_DIR)/.done_selected_gene_stats
$(RESPONSE_GENE_STATS_DONE): $(RESPONSE_SELECT_NETWORK_DONE)
	$(_start)
	$(_R) R/selected_gene_stats.r selected.gene.stats \
		ifn.genes=$(RESPONSE_NETWORK_GENES_SELECT) \
		poor.annotation.desc=$(UNIREF_POOR_ANNOTATION) \
		ifn.bg=$(UNIREF_STATS) \
		ofn=$(RESPONSE_NETWORK_GENES_SELECT_STATS)
	$(_end_touch)
selected_gene_stats: $(RESPONSE_GENE_STATS_DONE)

RESPONSE_GENE_WORDS_DONE?=$(EE_BASE_MAP_DIR)/.done_gene_words
$(RESPONSE_GENE_WORDS_DONE): $(RESPONSE_SELECT_NETWORK_DONE)
	$(_start)
	perl $(_md)/pl/gene_word.pl \
		$(RESPONSE_NETWORK_GENES_SELECT) \
		prot_desc \
		$(UNIREF_GENE_TAX_TABLE) \
		prot_desc \
		$(RESPONSE_GENE_WORD_TABLE) \
		$(RESPONSE_GENE_WORD_BACKTABLE) \
		$(RESPONSE_GENE_WORD_WHOLE) \
		$(UNIREF_POOR_ANNOTATION)
#	$(_end_touch)
gene_words: $(RESPONSE_GENE_WORDS_DONE)

# randomize network coords
rnetwork:
	rm -rf $(RESPONSE_COORD_NETWORK_DONE)
	@$(MAKE) network_coord

# use over response class
response_network: \
$(RESPONSE_CLUSTER_DONE) $(RESPONSE_INFO_DONE) $(RESPONSE_PATTERN_DONE) $(RESPONSE_GENES_DONE) $(SHORTLIST_ELEMENTS) \
$(RESPONSE_ANCHOR_ORDER_DONE) $(RESPONSE_ELEMENTS_FINAL_DONE) $(RESPONSE_COORD_NETWORK_DONE) $(RESPONSE_CLASSIFY_DONE) \
$(RESPONSE_GENE_STATS_DONE) $(MERGE_RESPONSE_GO_DONE)

##########################################################################################
# compare two networks
##########################################################################################

RESPONSE_NETWORK_COMPARE_DONE?=$(EE_COMPARE_MAPS_DIR)/.done_network_compare
$(RESPONSE_NETWORK_COMPARE_DONE):
	$(call _start,$(EE_COMPARE_MAPS_DIR))
	$(_R) R/network_compare.r network.compare \
		ifn.clusters=$(RESPONSE_CONTIG_CLUSTER) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.network1=$(RESPONSE_NETWORK1) \
		ifn.network2=$(RESPONSE_NETWORK2) \
		ifn.map1=$(RESPONSE_ANCHOR_MATRIX1) \
		ifn.map2=$(RESPONSE_ANCHOR_MATRIX2) \
		ifn.coverage1=$(COVERAGE_TABLE1) \
		ifn.coverage2=$(COVERAGE_TABLE2) \
		min.anchor.abundance=$(RESPONSE_NETWORK_MIN_ANCHOR_ABUNDNANCE) \
		min.element.abundance=$(RESPONSE_NETWORK_MIN_ELEMENT_ABUNDNANCE) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		min.enrichment=$(RESPONSE_NETWORK_MIN_ENRICHMENT) \
		separate.min.contacts=$(RESPONSE_NETWORK_SEPARATE_MIN_CONTACTS) \
		separate.max.enrichment=$(RESPONSE_NETWORK_SEPARATE_MAX_SCORE) \
		ofn.anchors=$(RESPONSE_ANCHOR_TABLE_COMPARE) \
		ofn.clusters=$(RESPONSE_CLUSTER_TABLE_COMPARE) \
		ofn.map=$(RESPONSE_MAP_COMPARE)
	$(_end_touch)
network_compare: $(RESPONSE_NETWORK_COMPARE_DONE)

RESPONSE_NETWORK_COMPARE_SELECT_DONE?=$(EE_COMPARE_MAPS_DIR)/.done_network_compare_select
$(RESPONSE_NETWORK_COMPARE_SELECT_DONE): $(RESPONSE_NETWORK_COMPARE_DONE)
	$(_start)
	$(_R) R/network_compare.r network.compare.select \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.map=$(RESPONSE_MAP_COMPARE) \
		ofn.map=$(RESPONSE_MAP_COMPARE_SELECT)
	$(_end_touch)
network_compare_select: $(RESPONSE_NETWORK_COMPARE_SELECT_DONE)

response_network_compare: $(RESPONSE_NETWORK_COMPARE_SELECT_DONE)

##########################################################################################
# deprecated
##########################################################################################

# plot per element the lca level of the hosts
RESPONSE_LEVEL_DONE?=$(EE_BASE_MAP_DIR)/.done_level_$(RESPONSE_NETWORK_LEVEL)
$(RESPONSE_LEVEL_DONE): $(RESPONSE_SELECT_NETWORK_DONE)
	$(_start)
	$(_R) R/element_host_analysis.r element.to.level \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.reps=$(SET_TAXA_REPS) \
		level=$(RESPONSE_NETWORK_LEVEL) \
		ofn.elements=$(RESPONSE_NETWORK_LEVEL_ELEMENTS) \
		ofn.matrix=$(RESPONSE_NETWORK_LEVEL_MATRIX) \
		ofn.summary=$(RESPONSE_NETWORK_LEVEL_SUMMARY)
	$(_end_touch)
element_level: $(RESPONSE_LEVEL_DONE)

RESPONSE_COORD_NETWORK_COMPARE_DONE?=$(EE_COMPARE_MAPS_DIR)/.done_network_coord_compare
$(RESPONSE_COORD_NETWORK_COMPARE_DONE):
	$(call _start,$(EE_COMPARE_MAPS_DIR))
	$(_R) R/network.r network.compare.coord \
		ifn1=$(RESPONSE_NETWORK1) \
		ifn2=$(RESPONSE_NETWORK2) \
		ofn=$(RESPONSE_NETWORK_COORDS_COMPARE)
	$(_end_touch)
response_network_compare_coords: $(RESPONSE_COORD_NETWORK_COMPARE_DONE)

# perl code converted to R code
#	perl $(_md)/pl/merge_ca_into_ea_matrix.pl \
#		$(CONTIG_TABLE) \
#		$(RESPONSE_CONTIG_CLUSTER) \
#		$(CA_MATRIX) \
#		$(EE_MATRIX) \
#		$(RESPONSE_ANCHOR_MATRIX)
