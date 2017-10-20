# contig coverage matrix
RESPONSE_CONTIG_DONE?=$(RESPONSE_DIR)/.done_contigs
$(RESPONSE_CONTIG_DONE):
	$(call _start,$(RESPONSE_DIR))
	$(_R) R/response.r compute.contig.response \
		ifn=$(CONTIG_TABLE) \
		min.detected=$(RESPONSE_CONTIG_OBSERVED_MIN_DETECTED) \
		assembly.dir=$(ASSEMBLY_DIR) \
		ids=$(TDATASETS) \
		map.tag=$(MAP_TAG) \
		ofn.observed=$(RESPONSE_CONTIG_OBSERVED) \
		ofn.expected=$(RESPONSE_CONTIG_EXPECTED) \
		ofn.norm=$(RESPONSE_CONTIG_NORM) \
		ofn.min.score=$(RESPONSE_CONTIG_NORM_DETECTION)
	$(_end_touch)
response_contig: $(RESPONSE_CONTIG_DONE)

################################################################################
# cluster contigs into elements by repsonse patterns
################################################################################

# contig clustering
RESPONSE_CLUSTER_DONE?=$(RESPONSE_CLUSTER_DIR)/.done_cluster
$(RESPONSE_CLUSTER_DONE): $(RESPONSE_CONTIG_DONE)
	$(call _start,$(RESPONSE_CLUSTER_DIR))
	$(_R) R/response_cluster.r cluster.contigs \
		ifn=$(RESPONSE_CONTIG_NORM) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		threshold=$(RESPONSE_CLUSTER_THRESHOLD) \
		ofn.table=$(RESPONSE_CONTIG_CLUSTER) \
		ofn.majors=$(RESPONSE_MAJORS)
	$(_end_touch)
response_cluster: $(RESPONSE_CLUSTER_DONE)

################################################################################
# various element properties
################################################################################

# cluster attributes
RESPONSE_INFO_DONE?=$(RESPONSE_CLUSTER_DIR)/.done_info
$(RESPONSE_INFO_DONE): $(RESPONSE_CLUSTER_DONE)
	$(_start)
	$(_R) R/response_info.r info \
		response.ifn=$(RESPONSE_CONTIG_NORM) \
		clusters.ifn=$(RESPONSE_CONTIG_CLUSTER) \
		contig.ifn=$(CONTIG_TABLE) \
		coverage.ifn=$(COVERAGE_TABLE) \
		gc.ifn=$(CONTIG_GC_TABLE) \
		gene.ifn=$(GENE_TABLE) \
		uniref.ifn=$(UNIREF_GENE_TAX_TABLE) \
		ca.ifn=$(CA_ANCHOR_CONTIGS) \
		ofn=$(RESPONSE_INFO)
	$(_end_touch)
response_info: $(RESPONSE_INFO_DONE)

# mean response
RESPONSE_PATTERN_DONE?=$(RESPONSE_CLUSTER_DIR)/.done_pattern
$(RESPONSE_PATTERN_DONE): $(RESPONSE_CLUSTER_DONE)
	$(call _start,$(RESPONSE_CLUSTER_DIR))
	$(_R) R/response.r pattern.response \
		ifn.norm=$(RESPONSE_CONTIG_NORM) \
		ifn.obs=$(RESPONSE_CONTIG_OBSERVED) \
		ifn.exp=$(RESPONSE_CONTIG_EXPECTED) \
		ifn.clusters=$(RESPONSE_CONTIG_CLUSTER) \
		ofn.obs=$(RESPONSE_PATTERN_OBS) \
		ofn.exp=$(RESPONSE_PATTERN_EXP) \
		ofn.mean=$(RESPONSE_PATTERN_MEAN) \
		ofn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ofn.top100=$(RESPONSE_PATTERN_TOP100) \
		ofn.bottom0=$(RESPONSE_PATTERN_BOTTOM0) \
		ofn.top95=$(RESPONSE_PATTERN_TOP95) \
		ofn.bottom05=$(RESPONSE_PATTERN_BOTTOM05) \
		ofn.top75=$(RESPONSE_PATTERN_TOP75) \
		ofn.bottom25=$(RESPONSE_PATTERN_BOTTOM25) \
		ofn.sd=$(RESPONSE_PATTERN_SD)
	$(_end_touch)
response_pattern: $(RESPONSE_PATTERN_DONE)

# genes per element
RESPONSE_GENES_DONE?=$(RESPONSE_CLUSTER_DIR)/.done_genes
$(RESPONSE_GENES_DONE): $(RESPONSE_CLUSTER_DONE)
	$(_start)
	$(_R) R/response.r element.genes \
		ifn.clusters=$(RESPONSE_CONTIG_CLUSTER) \
		ifn.genes=$(GENE_TABLE) \
		ifn.uniref=$(UNIREF_GENE_TAX_TABLE) \
		ofn=$(RESPONSE_GENES)
	$(_end_touch)
response_genes: $(RESPONSE_GENES_DONE)

RESPONSE_ANCHOR_ORDER_DONE?=$(RESPONSE_CLUSTER_DIR)/.done_host_cluster
$(RESPONSE_ANCHOR_ORDER_DONE): $(RESPONSE_CLUSTER_DONE)
	$(_start)
	$(_R) R/response.r response.anchor.order \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		type=$(RESPONSE_CLASS_TYPE) \
		class.count=$(RESPONSE_N_CLASSES) \
		max.height=$(RESPONSE_CLASS_MAX_HEIGHT) \
		base.ids=$(RESPONSE_BASE_IDS) \
		ofn=$(RESPONSE_ANCHOR_ORDER)
	$(_end_touch)
response_anchor_order: $(RESPONSE_ANCHOR_ORDER_DONE)

################################################################################
# classify genes
################################################################################

RESPONSE_CLASSIFY_DONE?=$(RESPONSE_CLUSTER_DIR)/.done_classify
$(RESPONSE_CLASSIFY_DONE): $(RESPONSE_CLUSTER_DONE)
	$(_start)
	$(_R) R/response_info.r classify \
		clusters.ifn=$(RESPONSE_CONTIG_CLUSTER) \
		contig.ifn=$(CONTIG_TABLE) \
		gene.ifn=$(GENE_TABLE) \
		uniref.ifn=$(UNIREF_GENE_TAX_TABLE) \
		ca.ifn=$(CA_ANCHOR_CONTIGS) \
		ofn=$(RESPONSE_CLASSIFY)
	$(_end_touch)
response_classify: $(RESPONSE_CLASSIFY_DONE)

####################################################################
# compile binaries
####################################################################

response_init: $(CORRELATE_CLUSTER_BINARY) $(COVERAGE_KMEANS_BINARY)
