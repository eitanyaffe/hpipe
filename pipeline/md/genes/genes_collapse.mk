############################################################################################################
# cluster genes by self blast
############################################################################################################

SET_FASTA_TYPE?=aa
ifeq ($(SET_FASTA_TYPE),aa)
SET_FASTA=$(SET_FASTA_AA)
BLAST_RULE=blast_aa
else
SET_FASTA=$(SET_FASTA_NT)
BLAST_RULE=ublast
endif

SET_CLUSTER_ID=$(CLUSTER_IDENTITY_THRESHOLD)_$(CLUSTER_COVERAGE_THRESHOLD)
SET_CLUSTER_BLAST?=$(SET_CLUSTER_DIR)/blast.result

SET_SELF_PARSE_DONE?=$(SET_CLUSTER_DIR)/.done_self_blast
$(SET_SELF_PARSE_DONE):
	$(_start)
	@$(MAKE) $(BLAST_RULE) \
		BLAST_DIR=$(SET_CLUSTER_DIR) \
		BLAST_QUERY_TABLE=$(SET_TABLE) \
		BLAST_QUERY_FASTA=$(SET_FASTA) \
		BLAST_TARGET_TABLE=$(SET_TABLE) \
		BLAST_TARGET_FASTA=$(SET_FASTA) \
		BLAST_RESULT=$(SET_CLUSTER_BLAST)
	$(_end_touch)

SET_CLUSTER_TABLE_DONE?=$(SET_CLUSTER_DIR)/.done_table_$(SET_CLUSTER_ID)
$(SET_CLUSTER_TABLE_DONE): $(SET_SELF_PARSE_DONE)
	$(_start)
	$(_md)/pl/cluster_genes_by_blast.pl \
		$(SET_TABLE) \
		$(SET_CLUSTER_BLAST) \
		$(CLUSTER_IDENTITY_THRESHOLD) \
		$(CLUSTER_COVERAGE_THRESHOLD) \
		$(SET_CLUSTER_PREFIX) \
		$(SET_CLUSTER_TABLE)
	$(_end_touch)

SET_CLUSTER_MAP_DONE?=$(SET_CLUSTER_DIR)/.done_map_$(SET_CLUSTER_ID)
$(SET_CLUSTER_MAP_DONE): $(SET_CLUSTER_TABLE_DONE)
	$(call _start,$(SET_CLUSTER_DIR))
	perl $(_md)/pl/collapse_map_single.pl \
		$(SET_CLUSTER_BLAST) \
		$(COLLAPSE_FIELD) \
		$(COLLAPSE_FUNCTION) \
		$(SET_CLUSTER_MAP)
	$(_end_touch)

cluster_sets: $(SET_CLUSTER_MAP_DONE)

