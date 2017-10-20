# split input fastq files
SPLIT_DONE=$(MAP_DIR)/.done_split
$(SPLIT_DONE):
	@echo '##' Assembly: $(ASSEMBLY_ID)
	@echo '##' Mapping library: $(MAP_LIB_ID)
	$(call _assert,MAP_LIB_ID MAP_INPUT PREPROC_FINAL_BASE)
	$(call _start,$(SPLIT_DIR))
	$(call _time,$(MAP_DIR),split) \
		$(_md)/pl/split_fastq.pl \
			$(SPLIT_DIR) \
			$(MAP_SPLIT_READS_PER_FILE) \
			$(MAP_SPLIT_TRIM) \
			$(MAP_SPLIT_READ_OFFSET) \
			$(MAP_READ_LENGTH) \
			$(MAP_INPUT_STAT) \
			$(MAP_INPUT)
	$(_end_touch)
split: $(SPLIT_DONE)

FILTER_QSUB_DIR?=$(QSUB_DATASET_DIR)/map_filter

# filter reads
FILTER_DONE?=$(MAP_DIR)/.done_filter
FILTER_SCRIPT?=$(_md)/pl/filter_map.pl
$(FILTER_DONE): $(VERIFY_PARSE_DONE)
	$(call _start,$(FILTER_DIR))
	mkdir -p $(FILTER_STAT_DIR)
	$(call _time,$(MAP_DIR),filter) $(_R) $(_md)/R/distrib_map.r distrib.filter \
		script=$(FILTER_SCRIPT) \
		idir=$(PARSE_DIR) \
		odir=$(FILTER_DIR) \
		sdir=$(FILTER_STAT_DIR) \
		sfile=$(FILTER_STAT_FILE) \
		min.score=$(MAP_MIN_QUALITY_SCORE) \
		min.length=$(MAP_MIN_LENGTH) \
		min.distance=$(MAP_MIN_EDIT_DISTANCE) \
		qsub.dir=$(FILTER_QSUB_DIR) \
		batch.max.jobs=$(NUM_MAP_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		jobname=mfilter
ifeq ($(REMOVE_TRANSIENT),T)
	echo cleaning transient files...
	rm -rf $(SPLIT_DIR) $(MAPPED_DIR) $(PARSE_DIR)
endif
	$(_end_touch)
map: $(FILTER_DONE)

PAIRED_QSUB_DIR?=$(QSUB_DATASET_DIR)/datasets/$(DATASET)/map_pair

# parse the bwa/sam output and create a paired coord table
PAIRED_DONE?=$(MAP_DIR)/.done_pair
PAIR_SCRIPT?=$(_md)/pl/pair_reads.pl
$(PAIRED_DONE): $(FILTER_DONE)
	$(call _start,$(PAIRED_DIR))
	mkdir -p $(PAIRED_STAT_DIR)
	$(call _time,$(MAP_DIR),pair) $(_R) $(_md)/R/distrib_map.r distrib.pair \
		script=$(PAIR_SCRIPT) \
		idir=$(FILTER_DIR) \
		odir=$(PAIRED_DIR) \
		sdir=$(PAIRED_STAT_DIR) \
		sfile=$(PAIRED_STAT_FILE) \
		qsub.dir=$(PAIRED_QSUB_DIR) \
		batch.max.jobs=$(NUM_MAP_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		jobname=mpair
	$(_end_touch)
pair: $(PAIRED_DONE)

# note we don't need map pairs here
COVERAGE_DONE?=$(COVERAGE_DIR)/.done
$(COVERAGE_DONE): $(FILTER_DONE)
	$(call _start,$(COVERAGE_DIR))
	$(call _time,$(MAP_DIR),coverage) $(_md)/pl/coverage.pl \
		$(CONTIG_TABLE) \
		$(FILTER_DIR) \
		$(MAP_BINSIZE) \
		$(COVERAGE_DIR)
	$(_end_touch)

# note we use all contigs when mapping but for coverage we use only the long contigs
$(COVERAGE_TABLE): $(COVERAGE_DONE)
	$(_start)
	$(_R) $(_md)/R/coverage_table.r coverage.table \
		table=$(CONTIG_TABLE) \
		idir=$(COVERAGE_DIR) \
		ofn=$@
	$(_end)
coverage: $(COVERAGE_TABLE)

###########################################################################
# plots
###########################################################################

make_plot_map_summary:
	$(_R) $(_md)/R/stats.r plot.stats.summary \
		adir=$(ASSEMBLY_DIR) \
		map.tag=$(MAP_TAG) \
		fdir=$(SET_FIGURE_DIR)/mapping \
		aid=$(ASSEMBLY_ID) \
		ids=$(LIB_IDS) \
		titles=$(LIB_TITLES) \
		filter.id=$(FILTER_ID)

.PHONY: coverage
map_all: pair coverage

###########################################################################
# plots (deprecated)
###########################################################################

make_plot_map:
	$(_R) $(_md)/R/stats.r plot.stats \
		input.ifn=$(MAP_INPUT_STAT) \
		filter.ifn=$(FILTER_STAT_FILE) \
		pair.ifn=$(PAIRED_STAT_FILE) \
		parse.ifn=$(PARSE_STAT_FILE) \
		id=$(LIB_ID) \
		fdir=$(MAP_FIGURE_DIR)/map_stats
