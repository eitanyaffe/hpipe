#####################################################################################################
# register module
#####################################################################################################

# preq variables
# LIB_ID: identifier of a paired end sequence library (Hi-C or regular shotgun)
# LIB_INPUT_DIRS: list of one or more directories with fastq files which contain one or more
#                 files matching *R[12].fastq
pr_vars:=OUTDIR SICKLE SEQPREP DECONSEQ_DIR

units:=preproc.mk

$(call _register_module,preproc,$(units),global,$(pr_vars))

#####################################################################################################
# lib pre-process
#####################################################################################################

# work here
PREPROC_DIR?=$(OUTDIR)/libs/$(LIB_ID)
LIB_COMPLEXITY_TABLE?=$(PREPROC_DIR)/complexity.table

# make final links here
PREPROC_FINAL_BASE?=$(OUTDIR)/libs_final
PREPROC_FINAL_DIR?=$(PREPROC_FINAL_BASE)/$(LIB_ID)_$(PREPROC_MODE)

# intermediate files
PREPROC_INTER_DIR?=$(TMPDIR)/libs/$(LIB_ID)

# sanger / illumina / solexa
SICKLE_TYPE?=sanger

# generate 1 or 2 output directories according to modes: simple / clean / both
#  simple: do not handle 3C ligations
#  clean: remove 3C all non-genomic 3C ligations
PREPROC_MODES?=simple

# simple / clean
PREPROC_MODE?=simple

# number of words per kmer for the clean_ligation utility tool
CL_KWORDS=1
CL_KSIZE=31
CL_SITE?=$(CUTTER_SITE)
CL_MIN_COVERAGE?=2
CL_MIN_RATIO?=2
CL_MIN_LENGTH?=50

# before seqprep we split the input files into files with up to SPLIT_SIZE sequences
SPLIT_SIZE?=10000000

# SeqPrep minimum length (each side must be over this threshold)
SEQPREP_LENGTH?=60

# set to -6 for phred+64 and to nothing for phred+33
SEQPREP_PHRED_FLAG?=-6

# sequence adapters
ADAPTOR1=AGATCGGAAGAGCAC
ADAPTOR2=AGATCGGAAGAGCGT

DECONSEQ_SCRIPT?=$(DECONSEQ_DIR)/deconseq.pl

# Alignment coverage threshold in percentage
DECONSEQ_COVERAGE?=10

# Alignment identity threshold in percentage
DECONSEQ_IDENTITY?=80

SICKLE_MAX_JOBS?=25
SEQPREP_MAX_JOBS?=25
DECONSEQ_MAX_JOBS?=50

# Name of deconseq database to use (human)
DECONSEQ_DBS?=hsref

