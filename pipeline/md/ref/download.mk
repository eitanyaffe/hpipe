
download_fasta?=$(GENOME_ORG_DIR)/.done
$(download_fasta):
	@rm -rf $(GENOME_ORG_DIR)
	$(call _start,$(GENOME_ORG_DIR))
	$(call _time,$(GENOME_DIR),download) $(_md)/pl/download_genomes.pl \
		$(DOWNLOAD_INPUT_TABLE) \
	 	$(GENEBANK_TABLE) \
	 	$(GENOME_ORG_DIR)
	$(_end_touch)

$(REF_ORIGINAL_FASTA): $(download_fasta)
	$(_start)
	cat `find $(GENOME_ORG_DIR) -name genome.fasta` > $@
	$(_end)

$(GENOME_CONTIG_TABLE): $(download_fasta)
	$(_start)
	$(_md)/pl/contig_summary.pl $(GENOME_ORG_DIR) genome.fasta $@
	$(_end)

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
