INDEX_PREFIX?=$(INDEX_DIR)/index
INDEX_DONE?=$(INDEX_DIR)/.done
$(INDEX_DONE):
	$(call _start,$(INDEX_DIR))
	$(BWA_BIN) index \
		-p $(INDEX_PREFIX) \
		$(MAP_SEQ_FILE)
	$(_end_touch)
index: $(INDEX_DONE)

MAP_DONE?=$(MAP_DIR)/.done_map_bwa
$(MAP_DONE): $(INDEX_DONE) $(SPLIT_DONE)
	$(call _start,$(MAPPED_DIR))
	$(if $(or $(_dry),$(wildcard $(SPLIT_DIR)/*)),,$(error no files found in $(SPLIT_DIR)))
	$(foreach IFN, $(wildcard $(SPLIT_DIR)/*), $(BWA_BIN) mem -T 1 \
		-t $(NTHREADS) \
		$(INDEX_PREFIX) \
		$(IFN) > \
		$(MAPPED_DIR)/$(notdir $(IFN)); $(ASSERT); )
	$(_end_touch)

PARSE_QSUB_DIR?=$(QSUB_DIR)/datasets/$(DATASET)/map_parse

# parse the bwa/sam output
PARSE_SCRIPT?=$(_md)/pl/parse_bwa_sam.pl
$(PARSE_DONE): $(MAP_DONE)
	$(call _start,$(PARSE_DIR))
	mkdir -p $(PARSE_STAT_DIR)
	$(call _time,$(PARSE_DIR)) \
		$(_R) $(_md)/R/distrib_parse_bwa.r distrib.parse.bwa \
		script=$(PARSE_SCRIPT) \
		idir=$(MAPPED_DIR) \
		odir=$(PARSE_DIR) \
		sdir=$(PARSE_STAT_DIR) \
		sfile=$(PARSE_STAT_FILE) \
		qsub.dir=$(PARSE_QSUB_DIR) \
		batch.max.jobs=$(NUM_MAP_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		jobname=mparse
	$(_end_touch)
parse: $(PARSE_DONE)

# verify parse using ref genome
SHOULD_VERIFY?=F
VERIFY_SCRIPT?=$(_md)/pl/verify_parse.pl
$(VERIFY_PARSE_DONE): $(PARSE_DONE)
	$(_start)
ifeq ($(SHOULD_VERIFY),T)
	$(call _time,$(PARSE_DIR),verify) \
		$(_R) $(_md)/R/distrib_parse_bwa.r distrib.verify.parse \
		script=$(VERIFY_SCRIPT) \
		ref=$(FULL_CONTIG_FILE) \
		idir=$(PARSE_DIR) \
		qsub.dir=$(PARSE_QSUB_DIR) \
		batch.max.jobs=$(NUM_MAP_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		jobname=mverify
endif
	$(_end_touch)

map_bwa: $(VERIFY_PARSE_DONE)
