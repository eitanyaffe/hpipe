PREDICT_DONE?=$(PREDICT_DIR)/.done
$(PREDICT_DONE):
	$(call _start,$(PREDICT_DIR))
	$(call _time,$(PREDICT_DIR)) $(MGM)/gmhmmp -r -m $(MGM)/MetaGeneMark_v1.mod \
		-A $(PREDICT_DIR)/aa_sequence.tmp \
		-D $(PREDICT_DIR)/nt_sequence.tmp \
		-o $(PREDICT_DIR)/table $(GENES_INPUT)
	$(_md)/pl/mgm_parse.pl \
		$(PREDICT_DIR)/table \
		pre \
		$(MGM_GENE_TABLE)
	cat $(PREDICT_DIR)/aa_sequence.tmp | \
		perl $(_md)/pl/mgm_process_fasta.pl pre $(MGM_GENE_TABLE) > \
		$(GENE_FASTA_AA)
	cat $(PREDICT_DIR)/nt_sequence.tmp | \
		perl $(_md)/pl/mgm_process_fasta.pl pre $(MGM_GENE_TABLE) > \
		$(GENE_FASTA_NT)
	$(_md)/pl/add_aa_length.pl \
		$(GENE_FASTA_AA) \
		$(MGM_GENE_TABLE) \
		$(GENE_TABLE)
	$(_end_touch)
predict_genes: $(PREDICT_DONE)

$(GENE_CLUSTER_TABLE): $(PREDICT_DONE)
	$(_start)
	@$(MAKE) cluster_sets \
		SET_CLUSTER_DIR=$(GENE_CLUSTER_DIR) \
		SET_TABLE=$(GENE_TABLE) \
		SET_FASTA_NT=$(GENE_FASTA_NT) \
		SET_FASTA_AA=$(GENE_FASTA_AA) \
		SET_CLUSTER_PREFIX=pg \
		SET_CLUSTER_TABLE=$(GENE_CLUSTER_TABLE) \
		SET_CLUSTER_MAP=$(GENE_CLUSTER_MAP)
	$(_end)
cluster_genes: $(GENE_CLUSTER_TABLE)

make_genes: $(PREDICT_DONE) $(GENE_CLUSTER_TABLE)
	@echo "DONE genes. CGENE_TABLE=$(GENE_CLUSTER_TABLE)"

.PHONY: cluster_genes make_genes
