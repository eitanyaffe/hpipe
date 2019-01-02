
$(call _register_module,global)

#####################################################################################################
# classes and paths
#####################################################################################################

# output directory
BASE_OUTDIR?=$(PWD)/output
OUTDIR?=$(BASE_OUTDIR)/$(PROJECT_ID)

# fast tmp directory used for intermediate lib preproc files
BASE_TMPDIR?=/tmp/$(shell whoami)
TMPDIR?=$(BASE_TMPDIR)/$(PROJECT_ID)

#####################################################################################################
# distributed run params
#####################################################################################################

QSUB_DIR?=$(OUTDIR)/distrib_shared
MAX_JOBS_FN?=$(CURDIR)/input/max_jobs
QSUB_DATASET_DIR?=$(QSUB_DIR)/assembly_$(ASSEMBLY_ID)/dataset_$(DATASET)
QSUB_ANCHOR_DIR?=$(QSUB_DATASET_DIR)/$(ANCHOR)_$(ANCHOR_LABEL)

# max number of threads
NTHREADS?=40

# default max number of jobs
MAX_JOBS?=40

# par: parallel, sge: sungrid
# [ currently only par checked ]
DTYPE?=par

#####################################################################################################
# paths need to be defined if working directly, without a docker image
#####################################################################################################

############## PREPROC ##############

# by default assume sickle is on path
SICKLE?=sickle

# SeqPrep tool (by default in path)
SEQPREP?=/usr/local/bin/SeqPrep

# human sequence removed by deconseq
DECONSEQ_DIR?=/relman01/shared/tools/deconseq-standalone-0.4.3/

############## ASSEMBLY ##############

# assembly
MEGAHIT?=/home/dethlefs/bin/megahit

# mummer path
MUMMER_DIR?=/home/eitany/work/download/MUMmer3.23

# bwa
BWA_BIN?=/home/eitany/work/download/bwa-0.7.12/bwa

# diamond
DIAMOND_BIN?=diamond

#####################################################################################################
# dataset entity
#####################################################################################################

# by default use the hi-c map
# LIB_ID?=$(HIC_LIB_ID)
# LIB_INPUT_DIRS?=$(HIC_READ_DIR)

# synonyms, kept for backward compatiblity
DATASET=$(LIB_ID)

DATASET_DIR?=$(ASSEMBLY_DIR)/datasets/$(DATASET)
DATASET_TMP_DIR?=$(ASSEMBLY_TMP_DIR)/datasets/$(DATASET)

# RE cutter site
CUTTER_SITE?=GATC
CUTTER_TITLE?=DpnII

#####################################################################################################
# anchor entity
#####################################################################################################

ANCHOR?=$(DATASET)

# change use to support multiple anchor inference parameter sets
ANCHOR_LABEL?=main
# ANCHOR_INFERENCE_DIR?=$(DATASET_DIR)/compute_anchors/$(ANCHOR_LABEL)

ANCHOR_DIR?=$(ASSEMBLY_DIR)/anchors/$(ANCHOR)_$(ANCHOR_LABEL)
ANCHOR_TABLE?=$(ANCHOR_DIR)/anchor_table
ORDER_DIR?=$(ANCHOR_DIR)/order
ANCHOR_CLUSTER_TABLE?=$(ORDER_DIR)/identity_$(SET_ID)

#####################################################################################################
# transient files
#####################################################################################################

# should remove transient files, keeping only result files. applies to:
# - preproc: keep final copy
# - mapping: keep filtered results
REMOVE_TRANSIENT?=T
