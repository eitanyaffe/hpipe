#####################################################################################################
# register module
#####################################################################################################

pr_vars:=MUMMER_DIR

units:=contig_matrix.mk contig_similarity.mk contig_cluster.mk qc.mk

$(call _register_module,ccluster,$(units),map,$(pr_vars))

#####################################################################################################
# variables
#####################################################################################################

### similarity ####

SIM_DIR?=$(ASSEMBLY_DIR)/similarity
SIM_RESULT_DIR?=$(SIM_DIR)/result

# mummer paths
MUMMER?=$(MUMMER_DIR)/mummer
NUCMER?=$(MUMMER_DIR)/nucmer
SHOWCOORD?=$(MUMMER_DIR)/show-coords

# intermediate dirs
SIM_SPLIT_DIR?=$(SIM_DIR)/split
SIM_QSUB_DIR?=$(QSUB_DIR)/assembly/$(ASSEMBLY_ID)/csim

# number of files
SIM_CHUNKS?=10

SIM_DONE?=$(SIM_DIR)/.done

### matrix ####

# extend similarity rects by this offset when masking out contacts between similar regions
CONTIG_MATRIX_SIMILARITY_OFFSET=2000

CONTIG_MATRIX_DIR?=$(DATASET_DIR)/contig_mat
CONTIG_CONTACTS?=$(CONTIG_MATRIX_DIR)/filtered_contacts
CONTIG_MASKED_CONTACTS?=$(CONTIG_MATRIX_DIR)/masked_contacts
CONTIG_MATRIX?=$(CONTIG_MATRIX_DIR)/table


CONTIG_MATRIX_FILTER_PARAMS?=$(CONTIG_MATRIX_DIR)/table_filter_params
CONTIG_MATRIX_FILTERED?=$(CONTIG_MATRIX_DIR)/table_filtered
CONTIG_MATRIX_STATS?=$(CONTIG_MATRIX_DIR)/stats

### cluster ####

# sample score if number over this threshold
CCLUSTER_MAX_CONTIG_PAIRS?=400

# cluster metric
#   contacts: log of number of contacts between contigs + 1 (log2(C_12+1)
#   shared: number of shared neighbours of the two contigs (S_12)
#   shared_corrected: shared neighbours corrected for expected shared function of marignals:
#       (1/N + S_12) / ((M1*M2)/N)
#   where N is total contigs and M1/2 is contig marginals
#   shared_corrected_reg: same as stSharedCorrected but regulated differently to prefer contigs with many share
#       (1 + S_12) / (1 + (M1*M2)/N)
CCLUSTER_SCORE_TYPE?=shared

# normalize score by cluster size
CCLUSTER_NORMALIZE_SIZE?=T

# stop below this score
CCLUSTER_MIN_SCORE?=50

# min marignal coverage of contig in terms of number of contig neighbors
CCLUSTER_MIN_COVERAGE?=20

# average percentage of shared neighbours over which we define a cell
CCLUSTER_CUTOFF?=90

# min elements in cluster
CCLUSTER_MIN_ELEMENTS?=20

# identifier of the initial anchor clustering
CLUSTER_ID?=score_$(CCLUSTER_MIN_SCORE)_marg_$(CCLUSTER_MIN_COVERAGE)_type_$(CCLUSTER_SCORE_TYPE)
CUTOFF_ID?=cutoff_$(CCLUSTER_CUTOFF)_minelements_$(CCLUSTER_MIN_ELEMENTS)

CCLUSTER_VERSION?=
CCLUSTER_BASE_DIR?=$(DATASET_DIR)/seed_anchors$(CCLUSTER_VERSION)
CCLUSTER_DIR?=$(CCLUSTER_BASE_DIR)/ccluster_$(CLUSTER_ID)
CCLUSTER_METRIC?=$(CCLUSTER_DIR)/metric
CCLUSTER_CONTIGS?=$(CCLUSTER_DIR)/contigs
CCLUSTER_TREE?=$(CCLUSTER_DIR)/tree
CCLUSTER_SCORES?=$(CCLUSTER_DIR)/scores

# sample metric in this ratio
CMETRIC_RATIO?=100

CELL_DIR?=$(CCLUSTER_BASE_DIR)/cassign_$(CUTOFF_ID)
INITIAL_ANCHOR_TABLE?=$(CELL_DIR)/contig_cell_table

# cluster/marked
CCLUSTER_COLOR_BY?=ccluster
CCLUSTER_FIGURE_DIR?=$(MAP_FIGURE_DIR)/ccluster


### QC ###

MAP_CONTACT_CLOSE_CIS_THRESHOLD=1000
MAP_CONTACT_FAR_CIS_THRESHOLD=2000
MAP_CONTACT_BINSIZE=1000
MAP_CONTACT_TABLE_MAX_READS?=-1
ifneq ($(MAP_CONTACT_TABLE_MAX_READS),-1)
MAP_CONTACT_TABLE?=$(DATASET_DIR)/ccontig_analysis_S$(MAP_CONTACT_TABLE_MAX_READS)
else
MAP_CONTACT_TABLE?=$(DATASET_DIR)/ccontig_analysis
endif

CIS_DECAY_BIN_LOG_START?=2
CIS_DECAY_BIN_LOG_END?=6
CIS_DECAY_BIN_LOG_STEP?=0.1
CIS_DECAY_GAP?=1000
CIS_DECAY_MAX_READS?=-1
ifneq ($(CIS_DECAY_MAX_READS),-1)
CIS_DECAY_DIR?=$(DATASET_DIR)/decay_S$(CIS_DECAY_MAX_READS)
else
CIS_DECAY_DIR?=$(DATASET_DIR)/decay
endif

# for QC purposes, define distal using this threshold (bp)
CIS_DECAY_DISTAL_THRESHOLD?=2000

CIS_DECAY_SUMMARY_YMAX?=0

CIS_DECAY_YLIM_N="0.000000001 1"

CIS_DECAY_TABLE?=$(CIS_DECAY_DIR)/cis_decay
CIS_DECAY_SUMMARY?=$(CIS_DECAY_DIR)/summary
