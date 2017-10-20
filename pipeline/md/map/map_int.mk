# DATASET: identifier of mapped library
pr_modules:=preproc assembly

units:=map.mk map_bwa.mk contacts.mk varisum.mk

$(call _register_module,map,$(units),$(pr_modules),)

#####################################################################################################
# mapping
#####################################################################################################

# currently support only bwa
MAP_TYPE?=bwa

# we use all the assembly as reference, including shorter discarded contigs
MAP_SEQ_FILE?=$(FULL_CONTIG_FILE)

# one index file for each assembly+mapper pair
INDEX_DIR?=$(ASSEMBLY_DIR)/map_index/$(MAP_TYPE)

# for mapping use by default only part of read
MAP_SPLIT_TRIM?=T
MAP_SPLIT_READ_OFFSET=10
MAP_READ_LENGTH=40
MAP_TAG?=$(MAP_SPLIT_TRIM)_$(MAP_SPLIT_READ_OFFSET)_$(MAP_READ_LENGTH)

# map uses two directories
MAP_DIR?=$(DATASET_DIR)/map_$(MAP_TAG)
MAP_TMPDIR?=$(DATASET_TMP_DIR)/map_$(MAP_TAG)

SPLIT_DIR?=$(MAP_DIR)/split
MAPPED_DIR?=$(MAP_DIR)/mapped_$(MAP_TYPE)
PARSE_DIR?=$(MAP_DIR)/parsed

# which lib to use for mapping
MAP_LIB_ID?=$(LIB_ID)

# input fastq files
MAP_INPUT?=$(wildcard $(PREPROC_FINAL_BASE)/$(MAP_LIB_ID)/R*)

# split input reads before mapping
MAP_SPLIT_READS_PER_FILE?=10000000

# number of parallel map jobs
NUM_MAP_JOBS?=20

# Phred of read
MAP_MIN_QUALITY_SCORE?=30

# in nt, the total length of all M segments in the CIGAR
MAP_MIN_LENGTH?=30

# The sam NM score as reported by bwa
MAP_MIN_EDIT_DISTANCE?=0

FILTER_ID?=s$(MAP_MIN_QUALITY_SCORE)_l$(MAP_MIN_LENGTH)_d$(MAP_MIN_EDIT_DISTANCE)

FILTER_DIR?=$(MAP_DIR)/filter_$(FILTER_ID)

PAIRED_DIR?=$(MAP_DIR)/pair_$(FILTER_ID)

PARSE_DONE?=$(MAP_DIR)/.done_parse_$(MAP_TYPE)
VERIFY_PARSE_DONE?=$(MAP_DIR)/.done_verify_$(MAP_TYPE)

#####################################################################################################
# stats
#####################################################################################################

# input reads
MAP_INPUT_STAT?=$(MAP_DIR)/stats_input

# parse reads
PARSE_STAT_DIR?=$(MAP_DIR)/parse_stat
PARSE_STAT_FILE?=$(MAP_DIR)/parse_stat.table

# filtering reads
FILTER_STAT_DIR?=$(MAP_DIR)/filter_stat_$(FILTER_ID)
FILTER_STAT_FILE?=$(MAP_DIR)/filter_stat_$(FILTER_ID).table

# pairing reads
PAIRED_STAT_DIR?=$(MAP_DIR)/pair_stat_$(FILTER_ID)
PAIRED_STAT_FILE?=$(MAP_DIR)/pair_stat_$(FILTER_ID).table

#####################################################################################################
# coverage
#####################################################################################################

MAP_BINSIZE?=100
COVERAGE_DIR?=$(DATASET_DIR)/coverage_$(MAP_TAG)/contigs
COVERAGE_TABLE=$(DATASET_DIR)/coverage_$(MAP_TAG)/table

#####################################################################################################
# genomic variation
#####################################################################################################

# margin from read sides, when inferring variants
VAR_MARGIN?=10
VAR_TABLE_NT?=$(MAP_DIR)/var_table_nt
VAR_TABLE_LINK?=$(MAP_DIR)/var_table_link

# variation summary
VAR_SUMMARY_PERCENTAGE_CUTOFF=1
VAR_SUMMARY_COUNT_CUTOFF=2
VAR_SUMMARY=$(MAP_DIR)/var_contig_summary
VAR_SUMMARY_BINS=$(MAP_DIR)/var_contig_summary_bins

# !!! kept in wrong directory - should be in response
VAR_SUMMARY_MULTI?=$(MAP_DIR)/contigs_variation_multi

#####################################################################################################
# contacts
#####################################################################################################

# cross coverage
CROSS_BINSIZE?=200
CROSS_MIN_DIST?=1000
CROSS_MAX_DIST?=4000
CROSS_COVERAGE_DIR?=$(DATASET_DIR)/cross_coverage/contigs_bin$(CROSS_BINSIZE)_min$(CROSS_MIN_DIST)_max$(CROSS_MAX_DIST)

INTRA_CONTIG_MIN_LENGTH?=10000
INTRA_FDIR?=$(FIGURE_DIR)/intra_contig
INTRA_CONTIG_CONTACT_DIR?=$(DATASET_DIR)/intra_contig_contacts

