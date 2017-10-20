
$(UNIREF_TABLE): $(GENE_REF_IFN)
	$(call _start,$(UNIREF_TABLE_DIR))
	$(_md)/pl/make_uniref_table.pl \
		$(GENE_REF_IFN) \
		$@
	$(_end)
uniref_table: $(UNIREF_TABLE)

$(UNIREF_TAX_LOOKUP):
	$(call _start,$(UNIREF_TABLE_DIR))
	$(call _time,$(UNIREF_TABLE_DIR),lookup_table) \
		$(_md)/pl/make_uniref_lookup.pl \
		$(GENE_REF_XML_IFN) \
		$@
	$(_end)
uniref_lookup: $(UNIREF_TAX_LOOKUP)

$(UNIREF_GENE_TABLE): $(UNIREF_TABLE) $(GENE_REF_IFN)
	cat $(GENE_REF_IFN) | $(_md)/pl/uniref_summary.pl > $@

UNIREF_BLAST_DONE?=$(UNIREF_DIR)/.done_blast
$(UNIREF_BLAST_DONE): $(UNIREF_GENE_TABLE)
	@mkdir -p $(UNIREF_DIAMOND_DB_DIR)
	$(call _start,$(UNIREF_DIR))
	@$(MAKE) blast_aa \
		BLAST_DIR=$(UNIREF_DIR) \
		BLAST_QUERY_TABLE=$(GENE_TABLE) \
		BLAST_QUERY_FASTA=$(GENE_FASTA_AA) \
		BLAST_TARGET_TABLE=$(UNIREF_GENE_TABLE) \
		BLAST_TARGET_FASTA=$(GENE_REF_IFN) \
		DIAMOND_INDEX=$(UNIREF_DIAMOND_DB) \
		BLAST_RESULT=$(UNIREF_RAW_OFN)
	$(_end_touch)
uniref_blast: $(UNIREF_BLAST_DONE)

UNIREF_UNIQ_DONE?=$(UNIREF_DIR)/.done_uniq
$(UNIREF_UNIQ_DONE): $(UNIREF_BLAST_DONE)
	$(_start)
	$(call _time,$(UNIREF_DIR),uniq) $(_md)/pl/blast_to_uniq.pl \
		$(UNIREF_RAW_OFN) \
		$(UNIREF_OFN_UNIQUE)
	$(_end_touch)
uniref_uniq: $(UNIREF_UNIQ_DONE)

$(UNIREF_TOP): $(UNIREF_BLAST_DONE)
	$(_start)
	perl $(_md)/pl/blast_to_all.pl \
		$(UNIREF_RAW_OFN) \
		$(TOP_IDENTITY_RATIO) \
		$(TOP_IDENTITY_DIFF) \
		$@
	$(_end)
uniref_top: $(UNIREF_TOP)

$(UNIREF_GENE_TAX_TABLE): $(UNIREF_UNIQ_DONE) $(UNIREF_TABLE) $(UNIREF_TAX_LOOKUP)
	$(_start)
	$(call _time,$(UNIREF_DIR),lookup_uniref) \
		$(_md)/pl/lookup_uniref.pl \
			$(UNIREF_TABLE) \
			$(UNIREF_OFN_UNIQUE) \
			$@
	$(_end)

UNIREF_STATS_DONE?=$(UNIREF_DIR)/.done_stats
$(UNIREF_STATS_DONE): $(UNIREF_GENE_TAX_TABLE)
	$(_start)
	$(_R) R/gene_stats.r gene.stats \
		ifn.genes=$(GENE_TABLE) \
		ifn.uniref=$(UNIREF_GENE_TAX_TABLE) \
		poor.annotation.desc=$(UNIREF_POOR_ANNOTATION) \
		ofn=$(UNIREF_STATS)
	$(_end_touch)
gene_stats: $(UNIREF_STATS_DONE)

################################################################################
# GO
################################################################################

GO_TREE_DONE?=$(GO_DIR)/.done_tree
$(GO_TREE_DONE):
	$(call _start,$(GO_DIR))
	$(_md)/pl/parse_go_obo.pl \
		$(GO_BASIC_OBO) \
		$(GO_TREE)
	$(_end_touch)
GO_tree: $(GO_TREE_DONE)

UNIPARC2UNIPROT_DONE?=$(GO_DIR)/.done_uniparc2uniprot
$(UNIPARC2UNIPROT_DONE):
	$(call _start,$(UNIREF_TABLE_DIR))
	$(call _time,$(UNIREF_TABLE_DIR),uniparc2uniprot_table) \
		$(_md)/pl/parse_uniparc_all.pl \
		$(UNIPARC_XML_IFN) \
		$(UNIPARC2UNIPROT_LOOKUP)
	$(_end_touch)
uniparc_lookup: $(UNIPARC2UNIPROT_DONE)

UNIPROT2GO_DONE?=$(GO_DIR)/.done_uniprot2go
$(UNIPROT2GO_DONE):
	$(call _start,$(UNIREF_TABLE_DIR))
	$(call _time,$(UNIREF_TABLE_DIR),uniprot_go_table) \
		$(_md)/pl/parse_goa_gaf.pl \
		$(GOA_UNIPROT_TABLE) \
		$(UNIPROT2GO_LOOKUP)
	$(_end_touch)
uniprot_go_table: $(UNIPROT2GO_DONE)

# GO annotation
UNIREF_GO_DONE?=$(UNIREF_DIR)/.done_GO
#$(UNIREF_GO_DONE): $(UNIREF_GENE_TAX_TABLE) $(UNIPROT_GO_TABLE)
$(UNIREF_GO_DONE): $(UNIPARC2UNIPROT_DONE) $(UNIPROT2GO_DONE)
	$(_start)
	$(_md)/pl/annotate_GO.pl \
		$(UNIREF_GENE_TAX_TABLE) \
		$(UNIREF_TAX_LOOKUP) \
		$(UNIPROT2GO_LOOKUP) \
		$(UNIPARC2UNIPROT_LOOKUP) \
		T \
		$(UNIREF_GENE_GO)
	$(_end_touch)
annotate_GO: $(UNIREF_GO_DONE)

################################################################################

make_uniref: $(UNIREF_GENE_TAX_TABLE) $(UNIREF_STATS_DONE) $(UNIREF_TOP) $(GO_TREE_DONE) $(UNIREF_GO_DONE)
