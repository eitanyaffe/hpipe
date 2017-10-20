############################################################################################################
# predict genes
############################################################################################################

GENOME_MGM_GENE_TABLE?=$(GENOME_PREDICT_DIR)/mgm_genes.txt
GENOME_BASIC_GENE_TABLE?=$(GENOME_PREDICT_DIR)/basic_genes.txt

GENOME_PREDICT_DONE?=$(GENOME_PREDICT_DIR)/.done
$(GENOME_PREDICT_DONE):
	$(call _start,$(GENOME_PREDICT_DIR))
	$(MGM)/gmhmmp -r -m $(MGM)/MetaGeneMark_v1.mod \
		-A $(GENOME_PREDICT_DIR)/proteins_aa.tmp \
		-D $(GENOME_PREDICT_DIR)/proteins_nt.tmp \
		-o $(GENOME_PREDICT_DIR)/mgm_table $(REF_FASTA)
	$(_md)/pl/mgm_parse.pl \
		$(GENOME_PREDICT_DIR)/mgm_table \
		ref \
		$(GENOME_MGM_GENE_TABLE)
	cat $(GENOME_PREDICT_DIR)/proteins_aa.tmp | \
		perl $(_md)/pl/mgm_process_fasta.pl ref $(GENOME_MGM_GENE_TABLE) > \
		$(GENOME_GENE_FASTA_AA)
	cat $(GENOME_PREDICT_DIR)/proteins_nt.tmp | \
		perl $(_md)/pl/mgm_process_fasta.pl ref $(GENOME_MGM_GENE_TABLE) > \
		$(GENOME_GENE_FASTA_NT)
	$(_md)/pl/add_aa_length.pl \
		$(GENOME_GENE_FASTA_AA) \
		$(GENOME_MGM_GENE_TABLE) \
		$(GENOME_BASIC_GENE_TABLE)
	$(_md)/pl/add_genome.pl \
		$(GENOME_BASIC_GENE_TABLE) \
		$(GENOME_GENE_TABLE)
	$(_end_touch)
ref_genes_local: $(GENOME_PREDICT_DONE)

############################################################################################################
# cluster genome genes
############################################################################################################

COLLAPSE_REF_GENE_DONE?=$(REF_GENE_CLUSTER_DIR)/.done_collapse
$(COLLAPSE_REF_GENE_DONE): $(GENOME_PREDICT_DONE)
	$(_start)
	@$(MAKE) cluster_sets \
		SET_CLUSTER_DIR=$(REF_GENE_CLUSTER_DIR) \
		SET_TABLE=$(GENOME_GENE_TABLE) \
		SET_FASTA_NT=$(GENOME_GENE_FASTA_NT) \
		SET_FASTA_AA=$(GENOME_GENE_FASTA_AA) \
		SET_CLUSTER_PREFIX=rg \
		SET_CLUSTER_TABLE=$(REF_GENE_CLUSTER_TABLE) \
		SET_CLUSTER_MAP=$(REF_GENE_CLUSTER_MAP)
	$(_end_touch)
ref_cluster: $(COLLAPSE_REF_GENE_DONE)

############################################################################################################
# map predicted to ref and back
############################################################################################################

REF_R2P_DONE?=$(REF_R2P_DIR)/.done_r2p
$(REF_R2P_DONE):
	$(_start)
	@$(MAKE) blast_aa \
		BLAST_DIR=$(REF_R2P_DIR) \
		BLAST_QUERY_TABLE=$(GENOME_GENE_TABLE) \
		BLAST_QUERY_FASTA=$(GENOME_GENE_FASTA_AA) \
		BLAST_TARGET_TABLE=$(GENE_TABLE) \
		BLAST_TARGET_FASTA=$(GENE_FASTA_AA) \
		BLAST_RESULT=$(REF_R2P_TABLE)
	$(_end_touch)
r2p: $(REF_R2P_DONE)

REF_P2R_DONE?=$(REF_P2R_DIR)/.done_p2r
$(REF_P2R_DONE):
	$(_start)
	@$(MAKE) blast_aa \
		BLAST_DIR=$(REF_P2R_DIR) \
		BLAST_QUERY_TABLE=$(GENE_TABLE) \
		BLAST_QUERY_FASTA=$(GENE_FASTA_AA) \
		BLAST_TARGET_TABLE=$(GENOME_GENE_TABLE) \
		BLAST_TARGET_FASTA=$(GENOME_GENE_FASTA_AA) \
		BLAST_RESULT=$(REF_P2R_TABLE)
	$(_end_touch)
p2r: $(REF_P2R_DONE)

$(REF_GENE_MAP): $(REF_P2R_DONE) $(REF_R2P_DONE) $(COLLAPSE_REF_GENE_DONE) $(GENE_CLUSTER_TABLE)
	$(call _start,$(GENE_CLUSTER_DIR))
	perl $(_md)/pl/collapse_map.pl \
		$(REF_R2P_TABLE) \
		$(REF_P2R_TABLE) \
		$(COLLAPSE_FIELD) \
		$(COLLAPSE_FUNCTION) \
		$@
	$(_end)
ref_map: $(REF_GENE_MAP)

# identify chimeric contigs
CLASSIFY_DONE?=$(REF_COMPARE_DIR)/.done_gene_classify
$(CLASSIFY_DONE): $(REF_GENE_MAP)
	$(_start)
	perl $(_md)/pl/classify_assembly.pl \
		$(FULL_CONTIG_TABLE) \
		$(GENOME_GENE_TABLE) \
		$(GENE_TABLE) \
		$(REF_GENE_MAP) \
		$(CLASSIFY_ID_THRESHOLD_STRONG) \
		$(CLASSIFY_ID_THRESHOLD_WEAK) \
		$(CLASSIFY_REF_GENE_TABLE) \
		$(CLASSIFY_GENE_TABLE) \
		$(CLASSIFY_CONTIG_TABLE)
	$(_end_touch)
classify: $(CLASSIFY_DONE)

plot_classify:
	$(_start)
	$(_R) R/ref_classify.r plot.classify \
		ref.genes.ifn=$(CLASSIFY_REF_GENE_TABLE) \
		genes.ifn=$(CLASSIFY_GENE_TABLE) \
		contigs.ifn=$(CLASSIFY_CONTIG_TABLE) \
		order.ifn=$(REF_CLUSTER_TABLE) \
		fdir=$(CLASSIFY_FDIR)
	$(_end)
# make all
make_ref: $(COLLAPSE_REF_GENE_DONE) $(REF_GENE_MAP) $(GENE_FATE_TABLE_DONE) $(CLASSIFY_DONE)
	@echo ref genome cluster: $(REF_GENE_CLUSTER_TABLE)
	@echo ref genome map: $(REF_GENE_MAP)

.PHONY: ref ref_map ref_cluster plot_classify


