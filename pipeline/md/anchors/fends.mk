#########################################################################
# fend table
#########################################################################

# make fragment table
FENDS_BASIC_DONE?=$(FENDS_ASSEMBLY_DIR)/.done_basic
$(FENDS_BASIC_DONE):
	$(call _start,$(FENDS_ASSEMBLY_DIR))
	$(_md)/pl/extract_sites.pl \
		$(CUTTER_SITE) \
		$(CONTIG_FILE) \
		$(FENDS_BASIC)
	$(_end_touch)
fends_basic: $(FENDS_BASIC_DONE)

# add contig coverage to fends
FENDS_COVERAGE_DONE?=$(FENDS_DIR)/.done_coverage
$(FENDS_COVERAGE_DONE): $(FENDS_BASIC_DONE)
	$(call _start,$(FENDS_DIR))
	$(_md)/pl/append_coverage.pl \
		$(COVERAGE_TABLE) \
		$(FENDS_MIN_RELATIVE_ABUNDANCE) \
		$(FENDS_BASIC) \
		$(FENDS_COVERAGE)
	$(_end_touch)
fends: $(FENDS_COVERAGE_DONE)

CONTIG_FENDS_DONE?=$(FENDS_DIR)/.done_contig_fends
$(CONTIG_FENDS_DONE): $(FENDS_BASIC_DONE)
	$(_start)
	$(_R) R/contig_fend_summary.r contig.fend.summary \
		ifn.fends=$(FENDS_BASIC) \
		ifn.contigs=$(CONTIG_TABLE) \
		ofn=$(CONTIG_FENDS_SUMMARY)
	$(_end_touch)
contig_fends: $(CONTIG_FENDS_DONE)

#########################################################################
# coords2fends
#########################################################################

COORD2FEND_QSUB_DIR?=$(QSUB_DIR)/datasets/$(DATASET)/coords2fends
MAT_DONE?=$(MAT_DIR)/.done

$(MAT_DONE): $(FENDS_COVERAGE_DONE)
	$(call _start,$(MAT_DIR))
	cp \
		$(FENDS_COVERAGE) \
		$(MAT_DIR)/fend.table
	$(call _time,$(MAT_DIR)) $(_md)/R/distrib_coord2fend.r \
		$(DATASET) \
	        $(CONTIG_CONTACTS) \
		$(MAT_DIR)/fend.table \
		$(MAP_READ_LENGTH) \
		$(MAP_SPLIT_READ_OFFSET) \
	        $(DISCARD_FACING_CIS_BELOW) \
		$(SEGMENT_LEN_THRESHOLD) \
		$(MAT_SPLIT_DIR) \
                $(MAT_DIR) \
		$(COORD2FEND_QSUB_DIR) \
		$(_md) \
		$(DTYPE) \
		$(MAX_JOBS_FN) \
		$(CONTIG_TABLE) \
		$(MAP_SPLIT_TRIM)
	$(_end_touch)
coords2fends: $(MAT_DONE)

contact_stats:
	mkdir -p $(ANCHOR_FIGURE_DIR)/stats
	cp $(MAT_DIR)/s0.fend.stats $(ANCHOR_FIGURE_DIR)/stats/s0_contact_stats.txt
