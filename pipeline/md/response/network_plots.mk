####################################################################
# Network plots
####################################################################

copy_network_files:
	mkdir -p $(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/files
	cp $(RESPONSE_NETWORK_SELECT) $(RESPONSE_NETWORK_ELEMENTS_SELECT) $(RESPONSE_NETWORK_GENES_SELECT) $(RESPONSE_NETWORK_GENES_SELECT_STATS) $(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/files

network_stats:
	$(_R) R/network_stats.r network.stats \
		ifn=$(RESPONSE_NETWORK_SELECT) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_stats

# the whole network
plot_network:
	$(_R) R/plot_network.r plot.network \
		lid=$(RESPONSE_LIB_ID) \
		ifn.reps=$(SET_TAXA_REPS) \
		ifn.majors=$(RESPONSE_MAJORS) \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.coords=$(RESPONSE_NETWORK_COORDS) \
		ifn.taxa.legend=$(SET_TAX_LEGEND) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.classify=$(RESPONSE_CLASSIFY) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network

# plot with new coordinates
replot_network:
	rm -rf $(RESPONSE_COORD_NETWORK_DONE)
	$(MAKE) network_coord plot_network

# degree analysis
plot_network_node_degree:
	$(_R) R/plot_network_props.r plot.node.degree \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.majors=$(RESPONSE_MAJORS) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_node_degree

plot_network_host2host_degree:
	$(_R) R/plot_network_props.r plot.host2host.degree \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.majors=$(RESPONSE_MAJORS) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_node_degree

plot_elements:
	$(_R) R/plot_elements.r plot.elements \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/elements

plot_network_chart:
	$(_R) R/plot_element_matrix.r plot.element.matrix \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.cluster.class=$(RESPONSE_CLASSIFY) \
		ifn.taxa.legend=$(SET_TAX_LEGEND) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_chart

plot_host_sharing:
	$(_R) R/plot_host_sharing.r plot.host.sharing \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/host_sharing

plot_element_matrix_summary:
	$(_R) R/plot_element_matrix_summary.r plot.element.matrix.summary \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.majors=$(RESPONSE_MAJORS) \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.cluster.class=$(RESPONSE_CLASSIFY) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/element_matrix_summary

plot_element_classify:
	$(_R) R/plot_element_class.r plot.class \
		ifn.class=$(RESPONSE_CLASSIFY) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/element_class

plot_diameter_analysis:
	$(_R) R/plot_element_matrix_summary.r plot.diameter.analysis \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		mat.ifn=$(ANCHOR_MATRIX_TABLE) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/diameter_analysis

# element coverage profiles
plot_element_response_simple:
	$(_R) R/plot_element_response_simple.r plot.element.response.simple \
		ifn.reps=$(SET_TAXA_REPS) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.norm=$(RESPONSE_CONTIG_NORM) \
		ifn.contigs=$(RESPONSE_CONTIG_CLUSTER) \
		ifn.bottom5=$(RESPONSE_PATTERN_BOTTOM05) \
		ifn.bottom25=$(RESPONSE_PATTERN_BOTTOM25) \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.top75=$(RESPONSE_PATTERN_TOP75) \
		ifn.top95=$(RESPONSE_PATTERN_TOP95) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		labels=$(TLABELS) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_element_response_simple

plot_element_response:
	$(_R) R/plot_element_response.r plot.elements \
		lid=$(RESPONSE_LIB_ID) \
		ifn.reps=$(SET_TAXA_REPS) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.map=$(RESPONSE_ANCHOR_MATRIX) \
		ifn.bottom0=$(RESPONSE_PATTERN_BOTTOM0) \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.top100=$(RESPONSE_PATTERN_TOP100) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		labels=$(TLABELS) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_element_response

plot_single_host_element_response:
	$(_R) R/single_host.r plot.single.host.element.detailed \
		lid=$(RESPONSE_LIB_ID) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS) \
		ifn.norm=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		base.min.correlation=$(RESPONSE_MIN_BASE_CORRELATION) \
		base.ids=$(RESPONSE_BASE_IDS) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		labels=$(TLABELS) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/single_host

plot_single_host_element_response_summary:
	$(_R) R/single_host.r plot.single.host.elements \
		lid=$(RESPONSE_LIB_ID) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS) \
		ifn.norm=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		labels=$(TLABELS) \
		ifn.classify=$(RESPONSE_CLASSIFY) \
		ifn.taxa.legend=$(SET_TAX_LEGEND) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/single_host_summary

plot_single_host_element_trend:
	$(_R) R/single_host_trend.r plot.single.host.element.trend \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS) \
		ifn.norm=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		base.min.correlation=$(RESPONSE_MIN_BASE_CORRELATION) \
		base.ids=$(RESPONSE_BASE_IDS) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		labels=$(TLABELS) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/single_host_trend

plot_network_element_support:
	$(_R) R/plot_network.r plot.network.element.support \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.map=$(RESPONSE_ANCHOR_MATRIX) \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.ca.matrix=$(CA_MATRIX) \
		ifn.clusters=$(RESPONSE_CONTIG_CLUSTER) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_element_support

plot_network_element_summary:
	$(_R) R/plot_network.r plot.network.element.summary \
		lid=$(RESPONSE_LIB_ID) \
		ifn.reps=$(SET_TAXA_REPS) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.obs=$(RESPONSE_PATTERN_OBS) \
		ifn.exp=$(RESPONSE_PATTERN_EXP) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		labels=$(TLABELS) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_element_summary

plot_hic_map:
	$(_R) R/plot_element_map.r plot.hic.map \
		lid=$(RESPONSE_LIB_ID) \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.map=$(RESPONSE_ANCHOR_MATRIX) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/hic_map

plot_gene_words:
	$(_R) R/plot_gene_annotation.r plot.words \
		ifn=$(RESPONSE_GENE_WORD_TABLE) \
		filter.words=$(RESPONSE_GENE_WORD_FILTER) \
		min.count=$(RESPONSE_GENE_WORD_MIN_COUNT) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/gene_diagram

# GO enrichments
plot_GO:
	$(_R) R/plot_gene_annotation.r plot.GO \
		ifn=$(RESPONSE_GO_MERGE) \
		min.count=$(RESPONSE_GO_MIN_COUNT) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/gene_diagram

# plot over response class

#response_plot_network_basic: network_stats copy_network_files plot_network plot_network_node_degree plot_network_host2host_degree plot_element_response plot_network_element_summary plot_network_chart plot_element_matrix_summary plot_diameter_analysis plot_host_sharing plot_hic_map

response_plot_network: plot_element_response_simple copy_network_files network_stats plot_single_host_element_response plot_network_element_support plot_network plot_network_node_degree plot_network_host2host_degree plot_element_response plot_network_element_summary plot_network_chart plot_element_matrix_summary plot_element_classify plot_diameter_analysis plot_single_host_element_trend plot_hic_map plot_host_sharing plot_gene_words plot_GO

####################################################################
# Network comparison
####################################################################

plot_network_compare_abundance:
	$(_R) R/plot_network_compare_analysis.r plot.abundance.summary \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.anchors=$(RESPONSE_ANCHOR_TABLE_COMPARE) \
		ifn.clusters=$(RESPONSE_CLUSTER_TABLE_COMPARE) \
		ifn.map=$(RESPONSE_MAP_COMPARE) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/abundance

plot_network_compare_analysis:
	$(_R) R/plot_network_compare_analysis.r plot.scatter.summary \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.anchors=$(RESPONSE_ANCHOR_TABLE_COMPARE) \
		ifn.clusters=$(RESPONSE_CLUSTER_TABLE_COMPARE) \
		ifn.map=$(RESPONSE_MAP_COMPARE) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/scatter_summary

plot_network_compare_details:
	$(_R) R/plot_network_compare_details.r plot.network.compare.details \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.anchors=$(RESPONSE_ANCHOR_TABLE_COMPARE) \
		ifn.clusters=$(RESPONSE_CLUSTER_TABLE_COMPARE) \
		ifn.map=$(RESPONSE_MAP_COMPARE) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		ifn.map.select=$(RESPONSE_MAP_COMPARE_SELECT) \
		legend1=$(RESPONSE_LEGEND1) \
		legend2=$(RESPONSE_LEGEND2) \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/details

plot_network_compare_matrix:
	$(_R) R/plot_network_compare_matrix.r plot.network.compare.matrix \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.anchors=$(RESPONSE_ANCHOR_TABLE_COMPARE) \
		ifn.clusters=$(RESPONSE_CLUSTER_TABLE_COMPARE) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		ifn.map=$(RESPONSE_MAP_COMPARE) \
		ifn.map.select=$(RESPONSE_MAP_COMPARE_SELECT) \
		zoom=F \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/matrix

plot_network_compare_matrix_zoom:
	$(_R) R/plot_network_compare_matrix.r plot.network.compare.matrix \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.anchors=$(RESPONSE_ANCHOR_TABLE_COMPARE) \
		ifn.clusters=$(RESPONSE_CLUSTER_TABLE_COMPARE) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		ifn.map=$(RESPONSE_MAP_COMPARE) \
		ifn.map.select=$(RESPONSE_MAP_COMPARE_SELECT) \
		zoom=T \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/matrix

network_compare_stats:
	$(_R) R/network_compare_stats.r network.stats \
		ifn.map.select=$(RESPONSE_MAP_COMPARE_SELECT) \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/stats

plot_network_compare_response:
	$(_R) R/plot_element_response.r plot.elements.compare \
		ifn.map=$(RESPONSE_MAP_COMPARE_SELECT) \
		lid=$(RESPONSE_LIB_ID) \
		ifn.reps=$(SET_TAXA_REPS) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.bottom0=$(RESPONSE_PATTERN_BOTTOM0) \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.top100=$(RESPONSE_PATTERN_TOP100) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		labels=$(TLABELS) \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/response

plot_network_compare: \
	plot_network_compare_abundance plot_network_compare_analysis plot_network_compare_details \
	plot_network_compare_matrix plot_network_compare_matrix_zoom \
	plot_network_compare_response network_compare_stats

####################################################################
# deprecated plots
####################################################################

# identity vs sharing
plot_identity_vs_sharing:
	$(_R) R/plot_network_props.r plot.identity.vs.sharing \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.majors=$(RESPONSE_MAJORS) \
		ifn.identity=$(ANCHOR_MATRIX_TABLE) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/identity_vs_sharing

plot_hic_map_compare_matrix_old:
	$(_R) R/plot_map_compare.r plot.hic.map.compare.matrix \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.network1=$(RESPONSE_NETWORK1) \
		ifn.network2=$(RESPONSE_NETWORK2) \
		ifn.map1=$(RESPONSE_ANCHOR_MATRIX1) \
		ifn.map2=$(RESPONSE_ANCHOR_MATRIX2) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		min.anchor.abundance=$(RESPONSE_NETWORK_MIN_ANCHOR_ABUNDNANCE) \
		min.element.abundance=$(RESPONSE_NETWORK_MIN_ELEMENT_ABUNDNANCE) \
		ifn.anchors.compare=$(RESPONSE_ANCHOR_TABLE_COMPARE) \
		ifn.clusters.compare=$(RESPONSE_CLUSTER_TABLE_COMPARE) \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/matrix

# plot per element the lca level of the hosts
NETWORK_REP_MIN_IDENTITY=90
NETWORK_REP_MIN_FRAC=50
plot_element_level:
	$(_R) R/element_host_analysis.r plot.element.host.lca \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.reps=$(SET_TAXA_REPS) \
		identity.threshold=$(NETWORK_REP_MIN_IDENTITY) \
		frac.threshold=$(NETWORK_REP_MIN_FRAC) \
		fdir=$(RES_FDIR)/element_level

plot_network_level:
	$(_R) R/plot_level.r plot.level.summary \
		ifn.elements=$(RESPONSE_NETWORK_LEVEL_ELEMENTS) \
		ifn.matrix=$(RESPONSE_NETWORK_LEVEL_MATRIX) \
		ifn.summary=$(RESPONSE_NETWORK_LEVEL_SUMMARY) \
		fdir=$(RES_FDIR)/datasets/$(RESPONSE_LIB_ID)/network_level_summary_$(RESPONSE_NETWORK_LEVEL)

plot_hic_map_compare_network:
	$(_R) R/plot_map_compare.r plot.hic.map.compare.network \
		ifn.coords=$(RESPONSE_NETWORK_COORDS_COMPARE) \
		ifn.response=$(RESPONSE_PATTERN_MEAN) \
		ifn.majors=$(RESPONSE_MAJORS) \
		ifn.network1=$(RESPONSE_NETWORK1) \
		ifn.network2=$(RESPONSE_NETWORK2) \
		ifn.map1=$(RESPONSE_ANCHOR_MATRIX1) \
		ifn.map2=$(RESPONSE_ANCHOR_MATRIX2) \
		min.contacts=$(RESPONSE_NETWORK_MIN_CONTACTS) \
		min.enrichment=$(RESPONSE_NETWORK_MIN_ENRICHMENT) \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)

plot_network_compare_chart:
	$(_R) R/plot_element_matrix.r plot.element.matrix.compare \
		ifn.network=$(RESPONSE_NETWORK_SELECT) \
		ifn.elements=$(RESPONSE_NETWORK_ELEMENTS_SELECT) \
		ifn.anchors=$(RESPONSE_ANCHOR_TABLE_COMPARE) \
		ifn.clusters=$(RESPONSE_CLUSTER_TABLE_COMPARE) \
		ifn.map=$(RESPONSE_MAP_COMPARE) \
		ifn.map.select=$(RESPONSE_MAP_COMPARE_SELECT) \
		ifn.cluster.class=$(RESPONSE_CLASSIFY) \
		ifn.taxa.legend=$(SET_TAX_LEGEND) \
		fdir=$(RES_FDIR)/datasets/hic_map_compare/$(RESPONSE_COMPARE_ID)/network_chart

