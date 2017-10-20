##########################################################################################
# get subtrees under
##########################################################################################

GR_TABLE_DONE?=$(GR_DIR)/.done_table
$(GR_TABLE_DONE):
	$(call _start,$(GR_DIR))
	$(_R) R/gr.r get.taxa.table \
		ifn.taxa=$(GR_TAXA) \
		ifn.tree=$(GR_TAXA_TREE) \
		ifn.genebank=$(GENEBANK_TABLE) \
		ofn=$(GR_TABLE)
	$(_end_touch)
gr_table: $(GR_TABLE_DONE)

GR_GOLD_TABLE_DONE?=$(GR_DIR)/.done_gold_table
$(GR_GOLD_TABLE_DONE):
	$(call _start,$(GR_DIR))
	$(_R) R/gr.r get.gold.table \
		ifn.gold.all=$(GR_GOLD_ALL) \
		ifn.gold.select=$(GR_GOLD_SELECT) \
		ifn.genebank=$(GENEBANK_TABLE) \
		ofn=$(GR_TABLE)
	$(_end_touch)
gr_gold_table: $(GR_GOLD_TABLE_DONE)

##########################################################################################
# sample a set of reference genomes
##########################################################################################

GR_SAMPLE_TABLE_DONE?=$(GR_SAMPLE_DIR)/.done_sample_table
$(GR_SAMPLE_TABLE_DONE): $(GR_GOLD_TABLE_DONE)
	$(call _start,$(GR_SAMPLE_DIR))
	$(_R) R/gr.r sample.taxa.table \
		ifn=$(GR_TABLE) \
		seed=$(GR_INDEX) \
		max.genomes=$(GR_MAX_GENOMES) \
		ofn=$(GR_SAMPLE_TABLE)
	$(_end_touch)
gr_sample_table: $(GR_SAMPLE_TABLE_DONE)

##########################################################################################
# download genomes
##########################################################################################

GR_DOWNLOAD_DONE?=$(GR_SAMPLE_DIR)/.done_download
$(GR_DOWNLOAD_DONE): $(GR_SAMPLE_TABLE_DONE)
	$(call _start,$(GR_DOWNLOAD_DIR))
	$(_md)/pl/download_genomes.pl \
		$(GR_SAMPLE_TABLE) \
	 	$(GENEBANK_TABLE) \
	 	$(GR_DOWNLOAD_DIR) \
		$(GR_DOWNLOAD_LOG)
	$(_end_touch)
gr_download: $(GR_DOWNLOAD_DONE)

$(GR_BASE_FASTA): $(GR_DOWNLOAD_DONE)
	$(_start)
	cat `find $(GR_DOWNLOAD_DIR) -name genome.fasta` > $@
	$(_end)
gr_base_fasta: $(GR_BASE_FASTA)

$(GR_GENOME_LOOKUP_TABLE): $(GR_DOWNLOAD_DONE)
	$(_start)
	$(_md)/pl/contig_summary.pl $(GR_DOWNLOAD_DIR) genome.fasta $@
	$(_end)
gr_contig_table: $(GR_GENOME_LOOKUP_TABLE)

GR_GENOME_DONE?=$(GR_SAMPLE_DIR)/.done_genome
$(GR_GENOME_DONE): $(GR_BASE_FASTA) $(GR_GENOME_LOOKUP_TABLE)
	$(_start)
	perl $(_md)/pl/generate_pseodo_genome.pl \
		$(GR_GENOME_LOOKUP_TABLE) \
		$(GR_BASE_FASTA) \
		$(GR_GENOME_TABLE) \
		$(GR_FASTA)
	$(_end_touch)
gr_genomes: $(GR_GENOME_DONE)

##########################################################################################
# predict genes
##########################################################################################

GR_GENES_DONE?=$(GR_SAMPLE_DIR)/.done_genes
$(GR_GENES_DONE): $(GR_GENOME_DONE)
	$(_start)
	$(MAKE) ref_genes_local GENOME_PREDICT_DIR=$(GR_GENE_DIR) REF_FASTA=$(GR_FASTA)
	$(_end_touch)
gr_genes: $(GR_GENES_DONE)

##########################################################################################
# cluster genes
##########################################################################################

# by nt
GR_CLUSTER_NT_DONE?=$(GR_SAMPLE_DIR)/.done_cluster_nt_$(GR_CLUSTER_ID_NT)
$(GR_CLUSTER_NT_DONE): $(GR_GENES_DONE)
	$(call _start,$(GR_GENE_CLUSTER_NT_DIR))
	@$(MAKE) cluster_sets \
		SET_CLUSTER_DIR=$(GR_GENE_BLAST_NT_DIR) \
		SET_TABLE=$(GR_GENE_TABLE) \
		SET_FASTA_NT=$(GR_GENE_FASTA_NT) \
		SET_FASTA_TYPE=nt \
		CLUSTER_IDENTITY_THRESHOLD=$(GR_CLUSTER_IDENTITY_THRESHOLD_NT) \
		CLUSTER_COVERAGE_THRESHOLD=$(GR_CLUSTER_COVERAGE_THRESHOLD_NT) \
		SET_CLUSTER_PREFIX=nt \
		SET_CLUSTER_TABLE=$(GR_GENE_CLUSTER_NT_TABLE) \
		SET_CLUSTER_MAP=$(GR_GENE_CLUSTER_NT_MAP)
	$(_end_touch)
gr_cluster_nt: $(GR_CLUSTER_NT_DONE)

# by aa
GR_CLUSTER_AA_DONE?=$(GR_SAMPLE_DIR)/.done_cluster_aa_$(GR_CLUSTER_ID_AA)
$(GR_CLUSTER_AA_DONE): $(GR_GENES_DONE)
	$(call _start,$(GR_GENE_CLUSTER_AA_DIR))
	@$(MAKE) cluster_sets \
		SET_CLUSTER_DIR=$(GR_GENE_BLAST_AA_DIR) \
		SET_TABLE=$(GR_GENE_TABLE) \
		SET_FASTA_AA=$(GR_GENE_FASTA_AA) \
		SET_FASTA_TYPE=aa \
		CLUSTER_IDENTITY_THRESHOLD=$(GR_CLUSTER_IDENTITY_THRESHOLD_AA) \
		CLUSTER_COVERAGE_THRESHOLD=$(GR_CLUSTER_COVERAGE_THRESHOLD_AA) \
		SET_CLUSTER_PREFIX=aa \
		SET_CLUSTER_TABLE=$(GR_GENE_CLUSTER_AA_TABLE) \
		SET_CLUSTER_MAP=$(GR_GENE_CLUSTER_AA_MAP)
	$(_end_touch)
gr_cluster_aa: $(GR_CLUSTER_AA_DONE)

gr_cluster_genes: $(GR_CLUSTER_NT_DONE) $(GR_CLUSTER_AA_DONE)

##########################################################################################
# cluster ref genomes (by aa)
##########################################################################################

GR_CLUSTER_GENOMES_DONE?=$(GR_SAMPLE_DIR)/.done_cluster_genomes
$(GR_CLUSTER_GENOMES_DONE): $(GR_CLUSTER_AA_DONE)
	$(_start)
	@$(MAKE) m=anchors make_set_clusters \
		SET_CLUSTER_DIR=$(GR_CLUSTER_GENOME_DIR) \
		SET_CLUSTER_INPUT_MATRIX=$(GR_GENE_TABLE) \
		SET_CLUSTER_INPUT_FIELD=genome \
		SET_CLUSTER_INPUT_MAP=$(GR_GENE_CLUSTER_AA_MAP) \
		SET_CLUSTER_TABLE=$(GR_GENOME_CLUSTER_TABLE) \
		SET_CLUSTER_PREFIX=r \
		SET_GENES_TABLE=$(GR_GENE_TABLE) \
		SET_CLUSTER_QSUB_DIR=$(GR_QSUB_MAP_DIR) \
		SET_SELF=T \
		SET_TITLE=reference
	$(_end_touch)
gr_cluster_genomes: $(GR_CLUSTER_GENOMES_DONE)

gr_all: $(GR_CLUSTER_NT_DONE) $(GR_CLUSTER_AA_DONE) $(GR_CLUSTER_GENOMES_DONE)

##########################################################################################
# shared distribution
##########################################################################################

GR_SHARED_DONE?=$(GR_SAMPLE_DIR)/.done_shared
$(GR_SHARED_DONE): $(GR_CLUSTER_GENOMES_DONE) $(GR_INTERFERENCE_DONE)
	$(_start)
	$(_R) R/gr.r shared.distrib \
		ifn.genes=$(GR_GENE_CLUSTER_NT_TABLE) \
		ifn.interference=$(GR_INTERFERENCE_TABLE) \
		ifn.genomes=$(GR_GENOME_CLUSTER_TABLE) \
		ifn.taxa=$(GR_SAMPLE_TABLE) \
		ifn.lookup=$(GR_GENOME_LOOKUP_TABLE) \
		ofn.table=$(GR_SHARED_TABLE) \
		ofn.distrib=$(GR_SHARED_DISTRIB)
	$(_end_touch)
gr_shared: $(GR_SHARED_DONE)

##########################################################################################
# gene interference (using nt)
##########################################################################################

GR_INTERFERENCE_DONE?=$(GR_SAMPLE_DIR)/.done_interference
$(GR_INTERFERENCE_DONE): $(GR_CLUSTER_NT_DONE)
	$(_start)
	perl $(_md)/pl/gene_interference.pl \
		$(GR_GENE_CLUSTER_NT_MAP) \
		$(GR_INTERFERENCE_MIN_IDENTITY) \
		$(GR_INTERFERENCE_MAX_IDENTITY) \
		$(GR_INTERFERENCE_MIN_COVERAGE) \
		$(GR_INTERFERENCE_TABLE)
	$(_end)
gr_interference: $(GR_INTERFERENCE_DONE)
