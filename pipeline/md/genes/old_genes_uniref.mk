############################################################################################################
# blast genes vs a database with usearch
############################################################################################################

UNIREF_P2R_DONE?=$(REF_P2R_DIR)/.done__p2r
$(REF_P2R_DONE):
	$(_start)
	@$(MAKE) blast \
		BLAST_DIR=$(REF_P2R_DIR) \
		BLAST_QUERY_TABLE=$(GENE_TABLE) \
		BLAST_QUERY_FASTA=$(GENE_FASTA_AA) \
		BLAST_TARGET_TABLE=$(GENOME_GENE_TABLE) \
		BLAST_TARGET_FASTA=$(GENOME_GENE_FASTA_AA) \
		BLAST_RESULT=$(REF_P2R_TABLE)
	$(_end_touch)
p2r: $(REF_P2R_DONE)

USEARCH_DB_DONE?=$(USEARCH_DB_DIR)/.done
$(USEARCH_DB_DONE):
	$(call _start,$(USEARCH_DB_DIR))
	$(call _time,$(USEARCH_DB_DIR)) \
		$(USEARCH_BIN) -makeudb_ublast $(GENE_REF_IFN) \
			-output $(USEARCH_DB_DIR)/table.udb
	$(_end_touch)
usearch_db: $(USEARCH_DB_DONE)

USEARCH_DONE?=$(USEARCH_ODIR)/.done
$(USEARCH_DONE): $(USEARCH_DB_DONE)
	$(call _start,$(USEARCH_ODIR))
	$(call _time,$(USEARCH_ODIR),blast) \
		$(USEARCH_BIN) -ublast $(USEARCH_QUERY) \
			-db $(USEARCH_DB_DIR)/table.udb \
	          	-threads $(NTHREADS) \
			-evalue $(USEARCH_EVALUE) \
			-blast6out $(USEARCH_OFN)
	$(_end_touch)
usearch: $(USEARCH_DONE)

$(USEARCH_OFN_UNIQUE): $(USEARCH_DONE)
	$(_start)
	$(_md)/pl/blast6out_to_uniq_v8_uniref.pl \
		$(USEARCH_OFN) \
		$@
	$(_end)
usearch_uniq: $(USEARCH_OFN_UNIQUE)

$(USEARCH_TOP): $(USEARCH_DONE)
	$(_start)
	perl $(_md)/pl/blast6out_to_top_uniref.pl \
		$(USEARCH_OFN) \
		$(TOP_IDENTITY_RATIO) \
		$(TOP_IDENTITY_DIFF) \
		$@
	$(_end)
usearch_top: $(USEARCH_TOP)

$(UNIREF_TABLE):
	$(call _start,$(UNIREF_TABLE_DIR))
	$(_md)/pl/make_uniref_table.pl \
		$(GENE_REF_IFN) \
		$@
	$(_end)
uniref_table: $(UNIREF_TABLE)

$(UNIREF_TAX_LOOKUP):
	$(call _start,$(UNIREF_TABLE_DIR))
	$(call _time,$(UNIREF_TABLE_DIR),tax_lookup_table) \
		$(_md)/pl/make_uniref_tax_lookup.pl \
		$(GENE_REF_XML_IFN) \
		$@
	$(_end)
uniref_tax_lookup: $(UNIREF_TAX_LOOKUP)

$(UNIREF_GENE_TAX_TABLE): $(USEARCH_OFN_UNIQUE) $(UNIREF_TABLE) $(UNIREF_TAX_LOOKUP)
	$(_start)
	$(call _time,$(USEARCH_ODIR),lookup_uniref) \
		$(_md)/pl/lookup_uniref.pl \
			$(UNIREF_TABLE) \
			$(USEARCH_OFN_UNIQUE) \
			$@
	$(_end)
make_uniref: $(UNIREF_GENE_TAX_TABLE) $(USEARCH_TOP)

# optional: add classes to gene table
$(UNIREF_GENE_TAX_TABLE_CLASS): $(UNIREF_GENE_TAX_TABLE)
	$(_start)
	$(call _time,$(USEARCH_DB_DIR),classify) \
		$(_md)/pl/uniref_classify.pl \
			$(UNIREF_GENE_TAX_TABLE) \
			$(UNIREF_CLASS_TABLE) \
			$@
	$(_end)

uniref_class: $(UNIREF_GENE_TAX_TABLE_CLASS)
	@echo "DONE gene/taxa table: $(UNIREF_GENE_TAX_TABLE_CLASS)"


############################################################################################################
# UNIPROT
# UNIPROT_XML?=/relman01/shared/databases/UniProt/uniprot_sprot.xml

# TBD:
# - make lookup from uniref to uniprot
# - make lookup uniprot to go
# - download go tree

############################################################################################################

