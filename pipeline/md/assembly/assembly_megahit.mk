ifeq ($(RARIFY_FOLD),0)
INPUT_CMD=cat $(ASSEMBLY_INPUT_PATTERN)
else
INPUT_CMD=cat $(RARIFY_READS)
endif

ASSEMBLY_WORK_DIR?=$(ASSEMBLY_DIR)/work
$(FULL_CONTIG_FILE): $(RARIFY_DONE)
	$(call _start,$(ASSEMBLY_DIR))
	@rm -rf $(ASSEMBLY_WORK_DIR)
	$(call _time,$(ASSEMBLY_DIR)) \
		$(MEGAHIT) \
			-m 0.8 \
			-o $(ASSEMBLY_WORK_DIR) \
			--k-min $(ASSEMBLY_MIN_KMER) \
			--k-max $(ASSEMBLY_MAX_KMER) \
			-t $(NTHREADS) \
			--input-cmd "$(INPUT_CMD)" \
			$(MEGA_HIT_PARAMS)
	cp $(ASSEMBLY_WORK_DIR)/final.contigs.fa $@
	$(_end)
$(FULL_CONTIG_TABLE): $(FULL_CONTIG_FILE)
	cat $(FULL_CONTIG_FILE) | $(_md)/pl/fasta_summary.pl > $@
basic_assembly: $(FULL_CONTIG_TABLE)

