
ifeq ($(MINIA_UNITIG),T)
ASSEMBLY_WORK_DIR?=$(ASSEMBLY_DIR)/work
else
ASSEMBLY_WORK_DIR?=$(ASSEMBLY_DIR)/work
endif

ifeq ($(MINIA_UNITIG),T)
MINIA_PARAMS:=-starter simple -traversal unitig
else
MINIA_PARAMS:=-starter best -traversal contig
endif

MLOG_FILE?=$(ASSEMBLY_WORK_DIR)/.log

MINIA_ASSEMBLY_DONE?=$(ASSEMBLY_WORK_DIR)/.done
$(MINIA_ASSEMBLY_DONE): $(ASSEMBLY_INPUT_FILE_TABLE) $(RARIFY_DONE)
	$(call _start,$(ASSEMBLY_WORK_DIR))
	cd $(ASSEMBLY_WORK_DIR) && $(call _time,$(ASSEMBLY_WORK_DIR)) \
		$(MINIA) \
			-kmer-size $(MINIA_KSIZE) \
			-abundance-min $(MINIA_MIN_COVERAGE) \
			-max-memory $(MINIA_MAX_MEMORY) \
			-in $(ASSEMBLY_INPUT_FILE_TABLE) \
			-nb-cores $(NTHREADS) \
			-out-dir $(ASSEMBLY_WORK_DIR) \
			-out $(ASSEMBLY_WORK_DIR)/final \
			$(MINIA_PARAMS) > $(MLOG_FILE)
	rm -rf $(ASSEMBLY_WORK_DIR)/trash*
	$(_end_touch)

$(FULL_CONTIG_FILE): $(MINIA_ASSEMBLY_DONE)
	cat $(ASSEMBLY_WORK_DIR)/final.contigs.fa | $(_md)/pl/simplify_minia_contigs.pl > $@

$(FULL_CONTIG_TABLE): $(FULL_CONTIG_FILE)
	cat $(FULL_CONTIG_FILE) | $(_md)/pl/fasta_summary.pl > $@
basic_assembly: $(FULL_CONTIG_TABLE)

