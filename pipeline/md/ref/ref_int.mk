
#####################################################################################################
# register module
#####################################################################################################

units:=download.mk simulate.mk
preqs_var:=GENEBANK_TABLE DOWNLOAD_INPUT_TABLE
$(call _register_module,ref,$(units),global,$(preqs_var))

#####################################################################################################
# define module variables
#####################################################################################################

# root dir
GENOME_DIR?=$(OUTDIR)/ref_genomes/$(GENOME_ID)

GENOME_ORG_DIR?=$(GENOME_DIR)/orgs

# maps accession and id, use with DOWNLOAD_INPUT_TABLE to determine species
GENOME_CONTIG_TABLE?=$(GENOME_DIR)/contig_table

# fasta broken down to contigs
REF_ORIGINAL_FASTA?=$(GENOME_DIR)/united

# pseudo genome fasta, Ns removed
REF_FASTA?=$(DATA_DIR)/genome.fasta
REF_GENOME_TABLE?=$(DATA_DIR)/genome.table

#####################################################################################################
# simulation parameters
#####################################################################################################

# simulate data
DATA_DIR?=$(OUTDIR)/simulated_reads/$(GENOME_ID)

# abundance weights are in 1X units
ABUNDANCE_MIN?=1
ABUNDANCE_MAX?=10

ABUNDANCE_MIN_HIC?=1
ABUNDANCE_MAX_HIC?=5

# assembly parameters
BASE_XCOVERAGE?=10
BASE_MOLECULE_LENGTH?=200
INSERT?=800
INSERT_SD?=200
READ_LENGTH?=250

# hic cis decay
HIC_READ_COUNT?=2000000
HIC_POWER_EXPONENT=-1
HIC_NOISE_PROB?=0.01
HIC_DISTAL_PROB?=0.5

# simulate reads
SIM_ASSEMBLY_ID?=x$(BASE_XCOVERAGE)_l$(READ_LENGTH)_i$(INSERT)_isd$(INSERT_SD)
SIM_HIC_ID?=c$(HIC_READ_COUNT)_p$(HIC_POWER_EXPONENT)_t$(HIC_TRANS_PROB)_d$(HIC_DISTAL_PROB)

# output fasta files: assembly
ADIR?=$(DATA_DIR)/assembly/$(SIM_ASSEMBLY_ID)
ASSEMBLY_READ_DIR?=$(ADIR)/reads

# output fasta files: Hi-C
HIC_DIR?=$(DATA_DIR)/hic/$(SIM_HIC_ID)
HIC_READ_DIR?=$(HIC_DIR)/hic_reads

