#####################################################################################################
# register module
#####################################################################################################

units:=response.mk network.mk response_plots.mk network_plots.mk GO.mk contact_local_analysis.mk
$(call _register_module,response,$(units),global,)

#####################################################################################################

RESPONSE_ID?=full_temporal
TDATASETS?=\
cipro_full_S1 cipro_full_S2 cipro_full_S3 cipro_full_S4 cipro_full_S5 cipro_full_S6 cipro_full_S7 cipro_full_S8 cipro_full_S9 cipro_full_S10 cipro_full_S11 cipro_full_S12 \
cipro_full_S13 cipro_full_S14 cipro_full_S15 cipro_full_S16  cipro_full_S17 cipro_full_S18 cipro_full_S19 cipro_full_S20 cipro_full_S21 cipro_full_S22 cipro_full_S23 cipro_full_S24
TLABELS?=S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20 S21 S22 S23 S24
TDAYS=-28 -15 -8 -2 -1 +1 +2 +3 +4 +5 +6 +7 +8 +10 +12 +14 +16 +18 +22 +25 +32 +47 +79 +80
TZOOM_DAYS?=-14 +21
TDISTURB_DAYS?=0 5

# compute network using the RESPONSE_LIB_ID library
RESPONSE_LIB_ID?=$(LIB_ID)

########################################################################
# normalized response table
########################################################################

RESPONSE_DIR?=$(ASSEMBLY_DIR)/response/$(RESPONSE_ID)
RESPONSE_CONTIG_OBSERVED?=$(RESPONSE_DIR)/contigs_observed
RESPONSE_CONTIG_EXPECTED?=$(RESPONSE_DIR)/contigs_expected

# minimal o/e is inferred given the following observed over the shortest contig
RESPONSE_CONTIG_OBSERVED_MIN_DETECTED?=1
RESPONSE_CONTIG_NORM?=$(RESPONSE_DIR)/contigs_norm
RESPONSE_CONTIG_NORM_DETECTION?=$(RESPONSE_DIR)/contigs_norm_detection

# weighted mean (by contig length) anchor response
RESPONSE_ANCHOR_GLOBAL?=$(RESPONSE_DIR)/global_anchor_response

########################################################################
# clustering contigs by response
########################################################################

RESPONSE_VERSION?=v8
RESPONSE_CLUSTER_THRESHOLD?=0.9
RESPONSE_TAG?=pearson_$(RESPONSE_CLUSTER_THRESHOLD)_$(ANCHOR)_$(ANCHOR_LABEL)_$(RESPONSE_VERSION)

RESPONSE_CLUSTER_DIR?=$(RESPONSE_DIR)/$(RESPONSE_TAG)

RESPONSE_CONTIG_CLUSTER_INIT?=$(RESPONSE_CLUSTER_DIR)/clusters_init

# how many response classes
RESPONSE_BASE_IDS?="cipro_full_S1 cipro_full_S2 cipro_full_S3 cipro_full_S4 cipro_full_S5"
RESPONSE_DISTURB_IDS?="cipro_full_S6 cipro_full_S7 cipro_full_S8 cipro_full_S9 cipro_full_S10"
RESPONSE_N_CLASSES?=15
RESPONSE_CLASS_MAX_HEIGHT?=0.15
RESPONSE_CLASS_TYPE=count
RESPONSE_CONTIG_ORDER=$(RESPONSE_CLUSTER_DIR)/order

# major element per anchor
RESPONSE_MAJORS?=$(RESPONSE_CLUSTER_DIR)/anchor_majors

RESPONSE_CONTIG_CLUSTER?=$(RESPONSE_CLUSTER_DIR)/clusters

# response info: various attributes of clusters
RESPONSE_INFO?=$(RESPONSE_CLUSTER_DIR)/cluster_info

# classify clusters (phage/plasmid etc)
RESPONSE_CLASSIFY?=$(RESPONSE_CLUSTER_DIR)/cluster_classify

# order patterns by hclust
RESPONSE_ANCHOR_ORDER?=$(RESPONSE_CLUSTER_DIR)/anchor_order

RES_FDIR?=$(FIGURE_DIR)/response/$(RESPONSE_ID)/$(RESPONSE_TAG)

RESPONSE_PATTERN_OBS?=$(RESPONSE_CLUSTER_DIR)/pattern_obs
RESPONSE_PATTERN_EXP?=$(RESPONSE_CLUSTER_DIR)/pattern_exp
RESPONSE_PATTERN_MEAN?=$(RESPONSE_CLUSTER_DIR)/pattern_mean
RESPONSE_PATTERN_MEDIAN?=$(RESPONSE_CLUSTER_DIR)/pattern_median
RESPONSE_PATTERN_TOP95?=$(RESPONSE_CLUSTER_DIR)/pattern_top95
RESPONSE_PATTERN_TOP75?=$(RESPONSE_CLUSTER_DIR)/pattern_top75
RESPONSE_PATTERN_TOP100?=$(RESPONSE_CLUSTER_DIR)/pattern_top100
RESPONSE_PATTERN_BOTTOM0?=$(RESPONSE_CLUSTER_DIR)/pattern_bottom0
RESPONSE_PATTERN_BOTTOM05?=$(RESPONSE_CLUSTER_DIR)/pattern_bottom05
RESPONSE_PATTERN_BOTTOM25?=$(RESPONSE_CLUSTER_DIR)/pattern_bottom25
RESPONSE_PATTERN_SD?=$(RESPONSE_CLUSTER_DIR)/pattern_sd

RESPONSE_GENES?=$(RESPONSE_CLUSTER_DIR)/genes

#############################################################################
# Network supported by the RESPONSE_LIB_ID Hi-C map
#############################################################################

# element-anchor association definition
# e: expected, o; observed, score: log10(0/e)
# connected: (score >= connected_min_score) && (o >= connected_min_o)
# separated: (score <= separated_max_score) && (e*10^separated_max_score >= separated_min_o)

# connection: minimal number of supporting contacts
RESPONSE_NETWORK_MIN_CONTACTS?=10

# connection: minimal score
RESPONSE_NETWORK_MIN_ENRICHMENT?=1

# separation: minimal number of expected contacts
RESPONSE_NETWORK_SEPARATE_MIN_CONTACTS?=0.5

# separation: maximal score
RESPONSE_NETWORK_SEPARATE_MAX_SCORE?=0.33

EE_BASE_MAP_DIR?=$(RESPONSE_CLUSTER_DIR)/ee_matrix/$(RESPONSE_LIB_ID)
EE_MAP_DIR=$(EE_BASE_MAP_DIR)/map

EE_MAP_IN_TABLE?=$(RESPONSE_CONTIG_CLUSTER)

# element-anchor association matrix
RESPONSE_ANCHOR_MATRIX?=$(EE_MAP_DIR)/matrix

RESPONSE_ANCHOR_ELEMENTS?=$(EE_BASE_MAP_DIR)/anchor_elements

# result network
RESPONSE_NETWORK?=$(EE_BASE_MAP_DIR)/network
RESPONSE_NETWORK_ELEMENTS?=$(EE_BASE_MAP_DIR)/network_elements
RESPONSE_NETWORK_GENES?=$(EE_BASE_MAP_DIR)/network_genes

#############################################################################
# single host analysis
#############################################################################

RESPONSE_MIN_BASE_CORRELATION?=0.98

#############################################################################
# selected network
#############################################################################

# minimal number of associated anchors
NETWORK_MIN_HOSTS?=2

# discard if is correlated with any anchor above this threshold
NETWORK_MAX_PEARSON?=0.99
#NETWORK_MAX_PEARSON?=1.00

# min/max number of genes
NETWORK_MIN_GENES?=2

# minimal number of uniref genes, which are not hypothetical or uncharacterized
NETWORK_MIN_UNIREF_GENES?=1

# keep element only if shared between any pair of distal anchors below identity threshold
#NETWORK_MAX_DIAMETER?=90
NETWORK_MAX_DIAMETER?=100

# result selected network
RESPONSE_NETWORK_SELECT?=$(EE_BASE_MAP_DIR)/select_network.txt
RESPONSE_NETWORK_ELEMENTS_SELECT?=$(EE_BASE_MAP_DIR)/select_network_elements.txt
RESPONSE_NETWORK_GENES_SELECT?=$(EE_BASE_MAP_DIR)/select_network_genes.txt

RESPONSE_NETWORK_GENES_SELECT_STATS?=$(EE_BASE_MAP_DIR)/select_network_genes_stats.txt

# consolidate anchors to genus level
RESPONSE_NETWORK_LEVEL?=family
RESPONSE_NETWORK_LEVEL_MATRIX?=$(EE_BASE_MAP_DIR)/level_matrix_$(RESPONSE_NETWORK_LEVEL)
RESPONSE_NETWORK_LEVEL_ELEMENTS?=$(EE_BASE_MAP_DIR)/level_elements_$(RESPONSE_NETWORK_LEVEL)
RESPONSE_NETWORK_LEVEL_SUMMARY?=$(EE_BASE_MAP_DIR)/level_summary_$(RESPONSE_NETWORK_LEVEL)

# for plotting, supported network layout types: graphopt/drl
RESPONSE_NETWORK_LAYOUT_TYPE?=graphopt
RESPONSE_NETWORK_COORDS?=$(EE_BASE_MAP_DIR)/network_coords_$(RESPONSE_NETWORK_LAYOUT_TYPE)

#############################################################################
# network comparison
#############################################################################

RESPONSE_COMPARE_ID?=pre_vs_post

RDATASET1?=pre_lib_hic_u_simple
RDATASET2?=post_lib_hic_u_simple

EE_COMPARE_MAPS_DIR?=$(RESPONSE_CLUSTER_DIR)/ee_matrix_compare/$(RESPONSE_COMPARE_ID)

RESPONSE_NETWORK1?=$(call reval,RESPONSE_NETWORK_SELECT,RESPONSE_LIB_ID=$(RDATASET1))
RESPONSE_NETWORK2?=$(call reval,RESPONSE_NETWORK_SELECT,RESPONSE_LIB_ID=$(RDATASET2))

RESPONSE_ANCHOR_MATRIX1?=$(call reval,RESPONSE_ANCHOR_MATRIX,RESPONSE_LIB_ID=$(RDATASET1))
RESPONSE_ANCHOR_MATRIX2?=$(call reval,RESPONSE_ANCHOR_MATRIX,RESPONSE_LIB_ID=$(RDATASET2))

COVERAGE_TABLE1=$(call reval,COVERAGE_TABLE,LIB_ID=$(RDATASET1))
COVERAGE_TABLE2=$(call reval,COVERAGE_TABLE,LIB_ID=$(RDATASET2))

RESPONSE_ANCHOR_TABLE_COMPARE?=$(EE_COMPARE_MAPS_DIR)/anchors_compare
RESPONSE_CLUSTER_TABLE_COMPARE?=$(EE_COMPARE_MAPS_DIR)/clusters_compare
RESPONSE_MAP_COMPARE?=$(EE_COMPARE_MAPS_DIR)/map_compare

# RESPONSE_NETWORK_COMPARE_TABLE?=$(EE_COMPARE_MAPS_DIR)/network_compare

RESPONSE_NETWORK_COORDS_COMPARE?=$(EE_COMPARE_MAPS_DIR)/network_coords

# disregard anchors/elements below this post abundance
RESPONSE_NETWORK_MIN_ANCHOR_ABUNDNANCE=-1
RESPONSE_NETWORK_MIN_ELEMENT_ABUNDNANCE=-1

# select significant anchor-element changes
RESPONSE_MAP_COMPARE_SELECT?=$(EE_COMPARE_MAPS_DIR)/map_compare_select

RESPONSE_LEGEND1?=pre
RESPONSE_LEGEND2?=post

#############################################################################
# compare Hi-C to standard shotgun abundance
#############################################################################

HIC_DATASET?=pre_lib_hic_u_simple
SG_DATASET?=cipro_full_S3
HIC_COVERAGE_TABLE=$(call reval,COVERAGE_TABLE,LIB_ID=$(HIC_DATASET))
SG_COVERAGE_TABLE=$(call reval,COVERAGE_TABLE,LIB_ID=$(SG_DATASET))

#############################################################################
# gene word table
#############################################################################

RESPONSE_GENE_WORD_MIN_COUNT?=3
RESPONSE_GENE_WORD_TABLE?=$(EE_BASE_MAP_DIR)/gene_word_table.txt
RESPONSE_GENE_WORD_BACKTABLE?=$(EE_BASE_MAP_DIR)/gene_word_backtable.txt
RESPONSE_GENE_WORD_FILTER?="Chromosome Partitioning Associated Proteins Dependent Predicted Modification Alpha C D 4 N Putative Like Core Site MULTISPECIES DNA RNA Protein Domain Subunit Major Cluster Related Region And Type Transcriptional Regulatory System Family I II S Specific Related_cluster Of Site-specific"
RESPONSE_GENE_WORD_WHOLE?="single-stranded helix-turn-helix DNA-binding transcriptional_regulator related_cluster ABC_transporter 50S_ribosomal_protein 30S_ribosomal_protein"

#############################################################################
# GO enrichments
#############################################################################

RESPONSE_GO_MIN_IDENTITY?=50
RESPONSE_GO_TABLE?=$(EE_BASE_MAP_DIR)/network_genes_GO
RESPONSE_GO_SUMMARY_PREFIX?=$(EE_BASE_MAP_DIR)/element_GO

RESPONSE_GO_TABLE_CTRL?=$(EE_BASE_MAP_DIR)/network_genes_GO_ctrl
RESPONSE_GO_SUMMARY_CTRL_PREFIX?=$(RESPONSE_GO_SUMMARY_PREFIX)_ctrl

RESPONSE_GO_MERGE?=$(EE_BASE_MAP_DIR)/GO_merge
RESPONSE_GO_MIN_COUNT?=3

#############################################################################
# local contact analysis
#############################################################################

LC_DIR?=$(EE_BASE_MAP_DIR)/contact_analysis
LC_READS?=$(LC_DIR)/reads
