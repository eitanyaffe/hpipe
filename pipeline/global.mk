
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
# external files
#####################################################################################################

# NCBI files, downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/
# First line (which has an #) was removed manually
# GENEBANK_DIR?=/relman01/shared/databases/NCBI/Genomes/Feb_2017
GENEBANK_TABLE?=$(GENEBANK_DIR)/assembly_summary_genbank.txt.2

#####################################################################################################
# exterior tools
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
MEGAHIT?=/home/eitany/work/download/megahit-master/megahit
MINIA?=/home/eitany/work/minia/bin/minia

# mummer path
MUMMER_DIR?=/home/eitany/work/download/MUMmer3.23

# bowtie
BOWTIE_BASE=/usr/local/bin
BOWTIE_BUILD_BIN?=$(BOWTIE_BASE)/bowtie-build
BOWTIE_BIN?=$(BOWTIE_BASE)/bowtie

# bwa
BWA_BIN?=/home/eitany/work/download/bwa-0.7.12/bwa

# diamond
DIAMOND_BIN?=diamond

############## CHECKM ##############

CHECKM?=~/work/python/local2/bin/checkm


############## BLAST ##############

USEARCH_BIN?=usearch8

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
# uniref
#####################################################################################################

GENE_REF_IFN?=$(GENE_REF_BASEDIR)/uniref100.fasta
GENE_REF_XML_IFN?=$(GENE_REF_BASEDIR)/uniref100.xmlx
GENE_REF_ID?=uniref100_2015_12

#####################################################################################################
# transient files
#####################################################################################################

# should remove transient files, keeping only result files. applies to:
# - preproc: keep final copy
# - mapping: keep filtered results
REMOVE_TRANSIENT?=T

#####################################################################################################
# figures
#####################################################################################################

BASE_FIGURE_DIR?=$(BASE_OUTDIR)/figures/$(PROJECT_ID)

# for various comparisons
SET_TITLE?=basic
SET_FIGURE_DIR?=$(BASE_FIGURE_DIR)/set_compare/$(SET_TITLE)

# for assembly
FIGURE_DIR?=$(BASE_FIGURE_DIR)/assembly_$(ASSEMBLY_ID)

LIB_FIGURE_DIR?=$(FIGURE_DIR)/lib_$(LIB_ID)
MAP_FIGURE_DIR=$(LIB_FIGURE_DIR)
GENE_FIGURE_DIR?=$(FIGURE_DIR)/genes
ANCHOR_FIGURE_DIR?=$(MAP_FIGURE_DIR)/anchor_$(ANCHOR)_$(ANCHOR_LABEL)

PP_COLS?="\#771155 \#AA4488 \#CC99BB \#114477 \#4477AA \#77AADD \#117777 \#44AAAA \#77CCCC \#117744 \#44AA77 \#88CCAA \#777711 \#AAAA44 \#DDDD77 \#774411 \#AA7744 \#DDAA77 \#771122 \#AA4455 \#DD7788"
LIB_COLS?="\#771155 \#AA4488 \#CC99BB \#114477 \#4477AA \#77AADD \#117777 \#44AAAA \#77CCCC \#117744 \#44AA77 \#88CCAA \#777711 \#AAAA44 \#DDDD77 \#774411 \#AA7744 \#DDAA77 \#771122 \#AA4455 \#DD7788"
ASSEMBLY_COLS?="\#771155 \#AA4488 \#CC99BB \#114477 \#4477AA \#77AADD \#117777 \#44AAAA \#77CCCC \#117744 \#44AA77 \#88CCAA \#777711 \#AAAA44 \#DDDD77 \#774411 \#AA7744 \#DDAA77 \#771122 \#AA4455 \#DD7788"
