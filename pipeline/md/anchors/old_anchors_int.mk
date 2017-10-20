#####################################################################################################
# register module
#####################################################################################################

units=fends.mk model.mk trim_anchors.mk ca_map.mk \
set_cluster.mk anchor_cluster.mk ref_cluster.mk \
set_compare.mk

$(call _register_module,anchors,$(units),ccluster)

#####################################################################################################
# Fends
#####################################################################################################

# RE cutter site
CUTTER_SITE?=GATC
CUTTER_TITLE?=DpnII

# per assembly
FENDS_ASSEMBLY_DIR?=$(ASSEMBLY_DIR)/fends
FENDS_BASIC?=$(FENDS_ASSEMBLY_DIR)/$(CUTTER_TITLE).basic

# per dataset
FENDS_DIR?=$(DATASET_DIR)/fends

# min compared to median coverage
FENDS_MIN_RELATIVE_ABUNDANCE?=0.01
FENDS_COVERAGE?=$(FENDS_DIR)/$(CUTTER_TITLE).coverage

MAT_DIR?=$(DATASET_DIR)/mat
MAT_SPLIT_DIR?=$(DATASET_DIR)/mat_split
MAT_FILE?=$(MAT_DIR)/s0.mat

# discard facing reads where the ends face each other and are near.
# These probably lack a ligation junction
DISCARD_FACING_CIS_BELOW?=2000

# discard reads for which the infered segment is short. Relevant for
# the 6 cutter Hi-C variant
SEGMENT_LEN_THRESHOLD?=1000

#####################################################################################################
# model parameters
#####################################################################################################

MDL_TYPE?=trans

# cutoff between far and close cis (bp)
CIS_THRESHOLD?=500000

MDL_FENDS?=$(FENDS_COVERAGE)
MDL_MAT?=$(MAT_FILE)

#####################################################################################################
# initial model (abundance only)
#####################################################################################################

# use the initial contig clustering
INIT_MDL_ANCHORS?=$(INITIAL_ANCHOR_TABLE)

# table with model parameters
INIT_MDL_MFN?=$(CURDIR)/input/models/abun.mdl

# output dir
INIT_MDL_DIR?=$(DATASET_DIR)/initial_model

# output prefix for all model related files
INIT_MODEL_PREFIX?=$(INIT_MDL_DIR)/$(MDL_TYPE)

INIT_MDL_DONE?=$(INIT_MDL_DIR)/.done

#####################################################################################################
# final model (abundance + fragment length)
#####################################################################################################

# use a final anchor table
FINAL_MDL_ANCHORS?=$(ANCHOR_TABLE)

# table with model parameters
FINAL_MDL_MFN?=$(CURDIR)/input/models/abun_len.mdl

# output dir
FINAL_MDL_DIR?=$(DATASET_DIR)/final_model

# output prefix for all model related files
FINAL_MODEL_PREFIX?=$(FINAL_MDL_DIR)/$(MDL_TYPE)

FINAL_MDL_DONE?=$(FINAL_MDL_DIR)/.done

#####################################################################################################
# set model params according to current stage
#####################################################################################################

MDL_STAGE?=initial

####################################
ifeq ($(MDL_STAGE),initial)
####################################

MDL_ANCHORS=$(INIT_MDL_ANCHORS)
MDL_ANCHOR_KEY=contig
MDL_ANCHOR_VALUE=cluster
MDL_MFN=$(INIT_MDL_MFN)

MDL_DIR=$(INIT_MDL_DIR)
MDL_QSUB=$(QSUB_DATASET_DIR)/initial_model

MDL_FDIR=$(MAP_FIGURE_DIR)/models/initial

MDL_ABUNDANCE_ONLY=T

BIAS_COLORS="black blue white red orange"
BIAS_BREAKS="-7 -2 0 2 7"

MDL_DONE=$(INIT_MDL_DONE)

####################################
else ifeq ($(MDL_STAGE),final)
####################################

MDL_ANCHORS=$(FINAL_MDL_ANCHORS)
MDL_ANCHOR_KEY=contig
MDL_ANCHOR_VALUE=anchor
MDL_MFN=$(FINAL_MDL_MFN)

MDL_DIR=$(FINAL_MDL_DIR)
MDL_QSUB=$(QSUB_DATASET_DIR)/final_model

MDL_FDIR=$(MAP_FIGURE_DIR)/models/final

MDL_ABUNDANCE_ONLY=F

BIAS_COLORS="black blue white red orange"
BIAS_BREAKS="-7 -2 0 2 7"

MDL_DONE=$(FINAL_MDL_DONE)

####################################
else
####################################

$(if $(MDL_STAGE),$(error unkown MDL_STAGE=$(MDL_STAGE)),$(error MDL_STAGE not defined))

####################################
endif
####################################


#####################################################################################################
# abundance model parameters
#####################################################################################################

# table with fends assigned to inclusive clusters
FENDS_CLUSTER_DIR?=$(FENDS_DIR)/$(CUTTER_TITLE)_$(CELL_ID)
FENDS_CLUSTER?=$(FENDS_CLUSTER_DIR)/fends

# number of abundance bins
ABUNDANCE_BINS?=100

ABUNDANCE_BINNED_FENDS?=$(INIT_MODEL_PREFIX).binned
ABUNDANCE_BIN_SIZE?=$(INIT_MODEL_PREFIX).total_counts
ABUNDANCE_COUNTS?=$(INIT_MODEL_PREFIX).counts

#####################################################################################################
# anchor trimming
#####################################################################################################

# paths
AMAP_OUT_DIR?=$(DATASET_DIR)/trim_anchors

# fend table and observed fend pairs
AMAP_IN_FENDS?=$(ABUNDANCE_BINNED_FENDS)
AMAP_IN_MAT?=$(MAT_FILE)

# model (for expected)
AMAP_IN_MODEL_PREFIX?=$(INIT_MODEL_PREFIX)

AMAP_IN_MODEL_FIELDS?=abundance

AMAP_NUM_SPLIT?=80

# only relevant for SGE
AMAP_REQ_MEM?=4000

# abundance only model file
AMAP_MFN?=$(INIT_MDL_MFN)

# minimal number of contacts
AMAP_MULTI_MIN_CONTACTS?=10

# number of bins when computing contig coverage
AMAP_MULTI_BIN_SIZE?=10

# min anchor size in bp
AMAP_MIN_LENGTH?=200000

# approximate FDR of contig-anchor trimming
AMAP_TRIM_FDR?=0.01

CELL_GENE_TAXA_TABLE?=$(AMAP_OUT_DIR)/cell_taxa


#####################################################################################################
# contig-anchor map (CA)
#####################################################################################################

# output dir
DATASET_ANCHOR_DIR=$(DATASET_DIR)/anchors/$(ANCHOR)
CA_MAP_DIR=$(DATASET_ANCHOR_DIR)/ca_matrix

# input
CA_MAP_IN_MODEL_PREFIX=$(FINAL_MODEL_PREFIX)
CA_MAP_IN_FENDS=$(CA_MAP_IN_MODEL_PREFIX).binned
CA_MAP_IN_MAT?=$(MAT_FILE)
CA_MAP_IN_MFN=$(FINAL_MDL_MFN)

CA_MAP_NUM_SPLIT?=80

# complete gene-anchor matrix
CA_MATRIX?=$(CA_MAP_DIR)/united

# to establish a contig-anchor association a pair must satify all conditions

# T = |anchors|*|contigs|
# N = |inter contacts|
# For each pair (c,a)_i assign a tail probability S_i assuming that obs_i is
# distributed as binomial B(N, P=exp_i/N), and require S_i < FDR / T
CA_ASSIGN_FDR?=0.000001

# min number of contacts required to consider contig-anchor or gene-anchor contact enrichment
CA_MIN_CONTACTS?=8

# minimal log10(obs/exp)
CA_MIN_ENRICHMENT?=1

# min number of anchor contigs required
CA_MIN_ANCHOR_CONTIGS?=4

# number of bins per contig, when computing contig coverage
CA_CONTIG_NBINS?=10

# min coverage of contig bins to associate all genes on contig to anchor
CA_CONTIG_COVERAGE?=0.5

# matrix limited to significant associations by thresholds set below
CA_ANCHOR_CONTIGS?=$(CA_MAP_DIR)/contigs
CA_ANCHOR_GENES?=$(CA_MAP_DIR)/genes

# table with various anchor info stats (contig median size, gc content etc)
ANCHOR_INFO_TABLE?=$(CA_MAP_DIR)/info_table
ANCHOR_INFO_TABLE_GENES?=$(CA_MAP_DIR)/info_table_genes
ANCHOR_INFO_TABLE_CONTIGS?=$(CA_MAP_DIR)/info_table_contigs

#####################################################################################################
# general gene-set comaprison parameters
#####################################################################################################

# cutoffs used when computing average identity
MEAN_IDENTITY_MIN_IDENTITY?=50
MEAN_IDENTITY_MIN_COVERAGE?=70

SET_ID?=$(COLLAPSE_BLAST_TYPE)_$(MEAN_IDENTITY_MIN_IDENTITY)_$(MEAN_IDENTITY_MIN_COVERAGE)

#####################################################################################################
# anchor cluster
#####################################################################################################

# working dir
ANCHOR_CLUSTER_DIR?=$(DATASET_ANCHOR_DIR)/cluster_identity_$(SET_ID)

# input set definition is CA
ANCHOR_CLUSTER_INPUT_MATRIX?=$(CA_ANCHOR_GENES)
ANCHOR_CLUSTER_INPUT_FIELD?=anchor

# input map
ANCHOR_CLUSTER_INPUT_MAP?=$(GENE_CLUSTER_MAP)

# output order file
ANCHOR_CLUSTER_TABLE?=$(ANCHOR_CLUSTER_DIR)/table.order

# qsub dir
ANCHOR_QSUB_MAP_DIR?=$(QSUB_DATASET_DIR)/cluster_anchor

# output figures here
ANCHOR_CLUSTER_FDIR?=$(ANCHOR_FIGURE_DIR)/cluster_identity_$(SET_ID)

#####################################################################################################
# ref cluster
#####################################################################################################

# working dir
REF_CLUSTER_DIR?=$(DATASET_ANCHOR_DIR)/ref_cluster_identity_$(SET_ID)

# input set
REF_CLUSTER_INPUT_MATRIX?=$(GENOME_GENE_TABLE)
REF_CLUSTER_INPUT_FIELD?=genome

# input map
REF_CLUSTER_INPUT_MAP?=$(REF_GENE_CLUSTER_MAP)

# output order file
REF_CLUSTER_TABLE?=$(REF_CLUSTER_DIR)/table.order

# qsub dir
REF_QSUB_MAP_DIR?=$(QSUB_DATASET_DIR)/cluster_ref

# output figures here
REF_CLUSTER_FDIR?=$(ANCHOR_FIGURE_DIR)/ref_cluster_identity_$(SET_ID)

#####################################################################################################
# comparison of two different gene sets (default: reference and anchor)
#####################################################################################################

SET_COMPARE_DIR?=$(DATASET_ANCHOR_DIR)/ref_compare_$(SET_ID)

# qsub dir
SET_COMPARE_QSUB?=$(QSUB_DATASET_DIR)/set_compare_ref_predicted

# gene set 1
REF_SET_INPUT1?=$(REF_CLUSTER_INPUT_MATRIX)
REF_SET_INPUT_GENES1?=$(REF_CLUSTER_INPUT_MATRIX)
REF_FIELD_INPUT1?=$(REF_CLUSTER_INPUT_FIELD)

# gene set 2
REF_SET_INPUT2?=$(ANCHOR_CLUSTER_INPUT_MATRIX)
REF_SET_INPUT_GENES2?=$(GENE_TABLE)
REF_FIELD_INPUT2?=$(ANCHOR_CLUSTER_INPUT_FIELD)

# map from 1 to 2
REF_MAP_INPUT?=$(REF_GENE_MAP)

REF_CLUSTER_TABLE1?=$(REF_CLUSTER_TABLE)
REF_CLUSTER_TABLE2?=$(ANCHOR_CLUSTER_TABLE)

REF_SET_TITLE1?=reference
REF_SET_TITLE2?=predicted

REF_SIMPLIFY_NAMES1?=F
REF_SIMPLIFY_NAMES2?=F

# mean/count
SET_COMPARE_METRIC?=count

SET_COMPARE_FDIR?=$(ANCHOR_FIGURE_DIR)/ref_compare_$(SET_ID)

# compare files
SET_COMPARE_MATRIX_DIR?=$(SET_COMPARE_DIR)/map.dir
SET_COMPARE_PROJECT_DIR?=$(SET_COMPARE_DIR)/project.dir
SET_COMPARE_MATRIX_TABLE?=$(SET_COMPARE_DIR)/summary.table

SET_MATCH_TABLE1?=$(SET_COMPARE_DIR)/match1
SET_MATCH_TABLE2?=$(SET_COMPARE_DIR)/match2

SET_FWD_TABLE?=$(SET_COMPARE_DIR)/fwd.table
SET_BCK_TABLE?=$(SET_COMPARE_DIR)/bck.table
SET_TMP_DIR?=$(SET_COMPARE_DIR)/tmp.dir

#####################################################################################################
# gene project
#####################################################################################################

GENE_PROJECT_DIR?=$(DATASET_ANCHOR_DIR)/projections

# discard weak hits for this map
GENE_PROJECT_MAP?=$(GENE_PROJECT_DIR)/map

# project genes forward from ref to predicted
GENE_PROJECT_FWD_SOURCE?=$(GENE_PROJECT_DIR)/gene_project_fwd_source
GENE_PROJECT_FWD_TARGET?=$(GENE_PROJECT_DIR)/gene_project_fwd_target

SET_FWD_SUMMARY_TABLE?=$(GENE_PROJECT_DIR)/fwd.summary
SET_BCK_SUMMARY_TABLE?=$(GENE_PROJECT_DIR)/bck.summary

SET_FWD_REPORT_TABLE?=$(GENE_PROJECT_DIR)/fwd.report
SET_BCK_REPORT_TABLE?=$(GENE_PROJECT_DIR)/bck.report

# project genes backwards from predicted to ref
GENE_PROJECT_BCK_SOURCE?=$(GENE_PROJECT_DIR)/gene_project_bck_source
GENE_PROJECT_BCK_TARGET?=$(GENE_PROJECT_DIR)/gene_project_bck_target

GENE_PROJECT_BACK?=$(GENE_PROJECT_DIR)/gene_project_back
GENE_PROJECT_IDENTITY_THRESHOLD=$(CLUSTER_IDENTITY_THRESHOLD)
GENE_PROJECT_COVERAGE_THRESHOLD=$(CLUSTER_COVERAGE_THRESHOLD)

GENE_PROJECT_FDIR?=$(ANCHOR_FIGURE_DIR)/projections
