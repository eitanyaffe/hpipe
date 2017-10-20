#####################################################################################################
# register module
#####################################################################################################

units=genes_predict.mk genes_uniref.mk genes_collapse.mk genes_ref.mk \
	genes_blast_aa.mk genes_blast_nt.mk genes_usearch.mk \
        genes_sample_ref.mk
var=GENE_REF_ID GENE_REF_IFN GENE_REF_XML_IFN
$(call _register_module,genes,$(units),anchors,$(var))

#####################################################################################################
# genes
#####################################################################################################

GENES_INPUT?=$(FULL_CONTIG_FILE)

# MetaGeneMarker
MGM?=/home/eitany/work/download/MetaGeneMark/MetaGeneMark_linux_64/mgm

# root dir for output
PREDICT_DIR?=$(ASSEMBLY_DIR)/genes

# temp table for mgm only
MGM_GENE_TABLE?=$(PREDICT_DIR)/mgm_genes.txt

# final gene table
GENE_TABLE?=$(PREDICT_DIR)/genes.txt
GENE_FASTA_AA?=$(PREDICT_DIR)/proteins_aa
GENE_FASTA_NT?=$(PREDICT_DIR)/proteins_nt

#####################################################################################################
# Gene Ontology
#####################################################################################################

# GO table from http://purl.obolibrary.org/obo/go/go-basic.obo
GO_BASIC_OBO?=/relman01/shared/databases/GO/go-basic.obo-2016-11-08

# UniProt to GO: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
GOA_UNIPROT_TABLE?=/relman01/shared/databases/GO/goa_uniprot/2016-10-29/goa_uniprot_all.gaf
UNIPROT2GO_LOOKUP?=/relman01/shared/databases/GO/goa_uniprot/2016-10-29/goa_uniprot_all.gaf.parsed

GO_ID?=2016_10
GO_DIR?=$(BASE_OUTDIR)/GO/$(GO_ID)
GO_TREE?=$(GO_DIR)/go_tree

#####################################################################################################
# uniref
#####################################################################################################

# UniRef from ftp://ftp.uniprot.org/pub/databases/uniprot/uniref
# GENE_REF_IFN?=/relman01/shared/databases/UniRef100/2014_11/uniref100.fasta
# GENE_REF_XML_IFN?=/relman01/shared/databases/UniRef100/2015_12/uniref100.xml
# GENE_REF_ID?=uniref100

# UniProt from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz
UNIPROT_XML_IFN?=/relman01/shared/databases/UniProt/versions/2016_10/uniprot_sprot.xml

# UniParc table: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_all.xml.gz
UNIPARC_XML_IFN?=/relman03/work/users/eitany/uniparc_all.xml

# uniref table
UNIREF_TABLE_DIR?=$(BASE_OUTDIR)/gene_catalog/$(GENE_REF_ID)
UNIREF_TABLE?=$(UNIREF_TABLE_DIR)/table
UNIREF_GENE_TABLE?=$(UNIREF_TABLE_DIR)/genes

# uniprot and taxa id lookup
UNIREF_TAX_LOOKUP?=$(UNIREF_TABLE_DIR)/tax_lookup

# uniparc to uniprot lookup
UNIPARC2UNIPROT_LOOKUP?=$(UNIREF_TABLE_DIR)/uniparc2uniprot_lookup

# diamond database
UNIREF_DIAMOND_DB_DIR?=$(BASE_OUTDIR)/diamond_db
UNIREF_DIAMOND_DB?=$(UNIREF_DIAMOND_DB_DIR)/$(GENE_REF_ID)

# uniref search result
UNIREF_ODIR_BASE?=$(ASSEMBLY_DIR)/genes_annotated
UNIREF_DIR?=$(UNIREF_ODIR_BASE)/$(GENE_REF_ID)
UNIREF_RAW_OFN?=$(UNIREF_DIR)/raw_table
UNIREF_OFN_UNIQUE?=$(UNIREF_DIR)/table_uniq
UNIREF_GENE_TAX_TABLE?=$(UNIREF_DIR)/table_uniq_taxa

# final annotated table, with GO
UNIREF_GENE_GO?=$(UNIREF_DIR)/table_GO

TOP_IDENTITY_RATIO=2
TOP_IDENTITY_DIFF=5
UNIREF_TOP?=$(UNIREF_DIR)/top_uniref_top

# gene general stats
UNIREF_POOR_ANNOTATION?="uncharacterized_protein hypothetical_protein MULTISPECIES:_hypothetical_protein Putative_uncharacterized_protein"
UNIREF_STATS?=$(UNIREF_DIR)/gene_stats

#####################################################################################################
# usearch (deprecated)
#####################################################################################################

# we blast using aa always
USEARCH_QUERY?=$(GENE_FASTA_AA)

USEARCH_EVALUE?=1e-5
USEARCH_ODIR_BASE?=$(ASSEMBLY_DIR)/genes_annotated
USEARCH_ODIR?=$(USEARCH_ODIR_BASE)/usearch_$(GENE_REF_ID)_$(USEARCH_EVALUE)
USEARCH_OFN?=$(USEARCH_ODIR)/raw_table
USEARCH_OFN_UNIQUE?=$(USEARCH_ODIR)/table_uniq

# index dir
USEARCH_DB_DIR?=$(BASE_OUTDIR)/usearch_index/$(GENE_REF_ID)

# with classification
UNIREF_CLASS_TABLE?=input/uniref_class
UNIREF_GENE_TAX_TABLE_CLASS?=$(USEARCH_ODIR)/table_uniq_taxa_class

#####################################################################################################
# cluster genes
#####################################################################################################

# general parameters
CLUSTER_IDENTITY_THRESHOLD?=95
CLUSTER_COVERAGE_THRESHOLD?=90

# if pair (gene1,gene2) appears more than once we select max|min over given field (coverage/evalue/identity)
COLLAPSE_FUNCTION?=max
COLLAPSE_FIELD?=identity
COLLAPSE_ID?=$(COLLAPSE_FIELD)_$(COLLAPSE_FUNCTION)

#####################################################################################################
# predicted genes
#####################################################################################################

GENE_CLUSTER_DIR?=$(ASSEMBLY_DIR)/gene_cluster
GENE_CLUSTER_TABLE?=$(GENE_CLUSTER_DIR)/cgene.table.$(CLUSTER_IDENTITY_THRESHOLD)_$(CLUSTER_COVERAGE_THRESHOLD)
GENE_CLUSTER_MAP?=$(GENE_CLUSTER_DIR)/gene_map_$(COLLAPSE_ID)

#####################################################################################################
# ref genes
#####################################################################################################

GENOME_PREDICT_DIR?=$(GENOME_DIR)/genes
GENOME_GENE_TABLE?=$(GENOME_PREDICT_DIR)/genes.txt
GENOME_GENE_FASTA_AA?=$(GENOME_PREDICT_DIR)/proteins_aa
GENOME_GENE_FASTA_NT?=$(GENOME_PREDICT_DIR)/proteins_nt

# ref gene clustering
REF_GENE_CLUSTER_DIR?=$(ASSEMBLY_DIR)/ref_gene_cluster
REF_GENE_CLUSTER_TABLE_INIT?=$(REF_GENE_CLUSTER_DIR)/init_cgene.table.$(CLUSTER_IDENTITY_THRESHOLD)_$(CLUSTER_COVERAGE_THRESHOLD)
REF_GENE_CLUSTER_TABLE?=$(REF_GENE_CLUSTER_DIR)/cgene.table.$(CLUSTER_IDENTITY_THRESHOLD)_$(CLUSTER_COVERAGE_THRESHOLD)
REF_GENE_CLUSTER_MAP?=$(REF_GENE_CLUSTER_DIR)/self_map_$(COLLAPSE_ID)

#####################################################################################################
# map from ref to predicted
#####################################################################################################

REF_COMPARE_DIR?=$(ASSEMBLY_DIR)/ref_compare

REF_GENE_MAP?=$(REF_COMPARE_DIR)/r2p_map_$(COLLAPSE_ID)

REF_R2P_DIR?=$(REF_COMPARE_DIR)/r2p
REF_R2P_TABLE?=$(REF_R2P_DIR)/map

REF_P2R_DIR?=$(REF_COMPARE_DIR)/p2r
REF_P2R_TABLE?=$(REF_P2R_DIR)/map

#####################################################################################################
# contig classify for chimeric identification
#####################################################################################################

CLASSIFY_ID_THRESHOLD_STRONG?=95
CLASSIFY_ID_THRESHOLD_WEAK?=70
CLASSIFY_REF_GENE_TABLE?=$(REF_COMPARE_DIR)/classify.ref_genes
CLASSIFY_GENE_TABLE?=$(REF_COMPARE_DIR)/classify.genes
CLASSIFY_CONTIG_TABLE?=$(REF_COMPARE_DIR)/classify.contigs
CLASSIFY_FDIR=$(FIGURE_DIR)/ref_classify

#####################################################################################################
# sample taxa space
#####################################################################################################

GR_VERSION?=v4
GR_LABEL?=basic
GR_ID?=$(GR_LABEL)_$(GR_VERSION)

GR_GOLD_ALL?=$(CURDIR)/input/GOLD/gold.tab
GR_GOLD_SELECT?=$(CURDIR)/input/GOLD/human_gut.tab

# taxa table, field: tax_id
GR_TAXA?=$(TAXA_SUBTREES_INPUT)

# taxa tree covering taxa table, fields: tax_id,parent_id
GR_TAXA_TREE?=$(TAXA_SUBTREES_TABLE)

# output dir
GR_DIR?=$(OUTDIR)/compare_ref_genomes/$(GR_ID)
GR_TABLE?=$(GR_DIR)/taxa_table

GR_INDEX?=1
GR_MAX_GENOMES?=500
GR_SAMPLE_ID?=i$(GR_INDEX)_n$(GR_MAX_GENOMES)
GR_SAMPLE_DIR?=$(GR_DIR)/sample_$(GR_SAMPLE_ID)
GR_SAMPLE_TABLE?=$(GR_SAMPLE_DIR)/taxa_table_sample

# download refs
GR_DOWNLOAD_DIR?=$(GR_SAMPLE_DIR)/genomes
GR_DOWNLOAD_LOG?=$(GR_SAMPLE_DIR)/download_failed_log
GR_BASE_FASTA?=$(GR_SAMPLE_DIR)/base.fasta
GR_GENOME_LOOKUP_TABLE?=$(GR_SAMPLE_DIR)/contig.table

# genomes
GR_GENOME_TABLE?=$(GR_SAMPLE_DIR)/genome.table
GR_FASTA?=$(GR_SAMPLE_DIR)/final.fasta

# genes
GR_GENE_DIR?=$(GR_SAMPLE_DIR)/genes
GR_GENE_TABLE?=$(GR_GENE_DIR)/genes.txt
GR_GENE_FASTA_AA?=$(GR_GENE_DIR)/proteins_aa
GR_GENE_FASTA_NT?=$(GR_GENE_DIR)/proteins_nt

# cluster genes by nt
GR_CLUSTER_IDENTITY_THRESHOLD_NT?=99
GR_CLUSTER_COVERAGE_THRESHOLD_NT?=100
GR_CLUSTER_ID_NT?=$(GR_CLUSTER_IDENTITY_THRESHOLD_NT)_$(GR_CLUSTER_COVERAGE_THRESHOLD_NT)
GR_GENE_BLAST_NT_DIR?=$(GR_SAMPLE_DIR)/blast_genes_nt
GR_GENE_CLUSTER_NT_DIR?=$(GR_SAMPLE_DIR)/cluster_genes_nt/$(GR_CLUSTER_ID_NT)
GR_GENE_CLUSTER_NT_TABLE?=$(GR_GENE_CLUSTER_NT_DIR)/genes.cluster
GR_GENE_CLUSTER_NT_MAP?=$(GR_GENE_CLUSTER_NT_DIR)/genes.map

# cluster genes by aa
GR_CLUSTER_IDENTITY_THRESHOLD_AA?=95
GR_CLUSTER_COVERAGE_THRESHOLD_AA?=90
GR_CLUSTER_ID_AA?=$(GR_CLUSTER_IDENTITY_THRESHOLD_AA)_$(GR_CLUSTER_COVERAGE_THRESHOLD_AA)
GR_GENE_BLAST_AA_DIR?=$(GR_SAMPLE_DIR)/blast_genes_aa
GR_GENE_CLUSTER_AA_DIR?=$(GR_SAMPLE_DIR)/cluster_genes_aa/$(GR_CLUSTER_ID_AA)
GR_GENE_CLUSTER_AA_TABLE?=$(GR_GENE_CLUSTER_AA_DIR)/genes.cluster
GR_GENE_CLUSTER_AA_MAP?=$(GR_GENE_CLUSTER_AA_DIR)/genes.map

# cluster ref genomes
GR_CLUSTER_GENOME_DIR?=$(GR_SAMPLE_DIR)/cluster_genomes
GR_GENOME_CLUSTER_TABLE?=$(GR_CLUSTER_GENOME_DIR)/table
GR_GENOME_SUMMARY_TABLE?=$(GR_CLUSTER_GENOME_DIR)/summary.table

GR_QSUB_MAP_DIR?=$(QSUB_DIR)/gr/$(GR_ID)/$(GR_SAMPLE_ID)

####################
# old
####################

# gene interference
GR_INTERFERENCE_MIN_IDENTITY?=98
GR_INTERFERENCE_MAX_IDENTITY?=99.5
GR_INTERFERENCE_MIN_COVERAGE?=90
GR_INTERFERENCE_TABLE?=$(GR_GENE_CLUSTER_NT_DIR)/gene_interference

# shared genes
GR_SHARED_TABLE?=$(GR_SAMPLE_DIR)/shared.distrib
GR_SHARED_DISTRIB?=$(GR_SAMPLE_DIR)/shared.table
