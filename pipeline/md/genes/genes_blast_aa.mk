
# params
DIAMOND_BLOCK_SIZE=20
DIAMOND_INDEX_CHUNKS=1
DIAMOND_THREADS=40
DIAMOND_EVALUE=0.001

# files
DIAMOND_INDEX?=$(BLAST_TARGET_FASTA).diamond_index
BLAST_RESULT_RAW?=$(BLAST_DIR)/blast.result.raw
BLAST_RESULT_SAM?=$(BLAST_DIR)/blast.result.sam

BLAST_INDEX_DONE?=$(DIAMOND_INDEX).done
$(BLAST_INDEX_DONE):
	$(call _start,$(BLAST_DIR))
	$(call _time,$(BLAST_DIR),index) $(DIAMOND_BIN) makedb \
		-c $(DIAMOND_INDEX_CHUNKS) \
		--in $(BLAST_TARGET_FASTA) \
		-p $(DIAMOND_THREADS) \
		-d $(DIAMOND_INDEX)
	$(_end_touch)

BLAST_DONE?=$(BLAST_DIR)/.done_raw
$(BLAST_DONE): $(BLAST_INDEX_DONE)
	$(call _start,$(BLAST_DIR))
	$(call _time,$(BLAST_DIR),blast) $(DIAMOND_BIN) blastp \
		-b $(DIAMOND_BLOCK_SIZE) \
		-c $(DIAMOND_INDEX_CHUNKS) \
		-d $(DIAMOND_INDEX) \
		-p $(DIAMOND_THREADS) \
		-q $(BLAST_QUERY_FASTA) \
		-e $(DIAMOND_EVALUE) \
		--sensitive \
		-a $(BLAST_RESULT_RAW)
	$(_end_touch)

BLAST_SAM_DONE?=$(BLAST_DIR)/.done_sam
$(BLAST_SAM_DONE): $(BLAST_DONE)
	$(_start)
	$(call _time,$(BLAST_DIR),view) $(DIAMOND_BIN) view \
		-a $(BLAST_RESULT_RAW) \
		-o $(BLAST_RESULT_SAM)
	$(_end_touch)

BLAST_PARSE_DONE?=$(BLAST_DIR)/.done_parse
$(BLAST_PARSE_DONE): $(BLAST_SAM_DONE)
	$(_start)
	$(call _time,$(BLAST_DIR),parse_sam) perl $(_md)/pl/sam_parse.pl \
		$(BLAST_RESULT_SAM) \
		$(BLAST_QUERY_TABLE) \
		$(BLAST_TARGET_TABLE) \
		aa \
		$(BLAST_RESULT)
	$(_end_touch)
blast_aa: $(BLAST_PARSE_DONE)
