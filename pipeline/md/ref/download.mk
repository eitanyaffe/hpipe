
DOWNLOAD_FASTA_DONE?=$(GENOME_ORG_DIR)/.done
$(DOWNLOAD_FASTA_DONE):
	@rm -rf $(GENOME_ORG_DIR)
	$(call _start,$(GENOME_ORG_DIR))
	$(call _assert,GENEBANK_DIR)
	$(call _time,$(GENOME_DIR),download) $(_md)/pl/download_genomes.pl \
		$(DOWNLOAD_INPUT_TABLE) \
	 	$(GENEBANK_TABLE) \
	 	$(GENOME_ORG_DIR)
	$(_end_touch)
download_base: $(DOWNLOAD_FASTA_DONE)

DOWNLOAD_INPUT_DONE?=$(GENOME_ORG_DIR)/.done_input
$(DOWNLOAD_INPUT_DONE): $(DOWNLOAD_FASTA_DONE)
	$(call _start,$(GENOME_INPUT_DIR))
	$(_md)/pl/download_process.pl \
		$(GENOME_ORG_DIR) \
		genome.fasta \
		$(TRUNCATE_DOWNLOADED_GENOMES) \
		$(TRUNCATE_DOWNLOADED_GENOMES_LINES) \
		$(GENOME_INPUT_DIR)
	$(_end_touch)
download_input: $(DOWNLOAD_INPUT_DONE)

$(REF_ORIGINAL_FASTA): $(DOWNLOAD_INPUT_DONE)
	$(_start)
	cat `find $(GENOME_INPUT_DIR) -type f` > $@
	$(_end)
download_merge: $(REF_ORIGINAL_FASTA)

$(GENOME_CONTIG_TABLE): $(DOWNLOAD_INPUT_DONE)
	$(_start)
	$(_md)/pl/contig_summary.pl $(GENOME_INPUT_DIR) $@
	$(_end)
download_summary: $(GENOME_CONTIG_TABLE)

##################################################################################################
# construct psuedo-genomes (concatenate all contigs and discard Ns)
##################################################################################################

REF_GENOME_TABLE_INIT?=$(DATA_DIR)/genome_init.table

REF_GENOME_DONE?=$(DATA_DIR)/.done_genome
$(REF_GENOME_DONE): $(GENOME_CONTIG_TABLE)
	$(call _start,$(DATA_DIR))
	perl $(_md)/pl/generate_pseodo_genome.pl \
		$(GENOME_CONTIG_TABLE) \
		$(REF_ORIGINAL_FASTA) \
		$(REF_GENOME_TABLE_INIT) \
		$(REF_FASTA)
	$(_end_touch)

$(REF_GENOME_TABLE): $(REF_GENOME_DONE)
	$(_start)
	$(_R) $(_md)/R/ref.r generate.abundance.table \
		ifn=$(REF_GENOME_TABLE_INIT) \
		abundance.min=$(ABUNDANCE_MIN) \
		abundance.max=$(ABUNDANCE_MAX) \
		abundance.min.hic=$(ABUNDANCE_MIN_HIC) \
		abundance.max.hic=$(ABUNDANCE_MAX_HIC) \
		ofn=$@
	$(_end)

############################################################################################################
# all
############################################################################################################

download: $(REF_ORIGINAL_FASTA) $(REF_GENOME_TABLE)
