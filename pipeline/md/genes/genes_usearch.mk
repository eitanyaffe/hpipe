USEARCH_EVALUE=0.001
USEARCH_THREADS?=40
USEARCH_TYPE?=nt
USEARCH_ID?=0.9

USEARCH_RESULT_SAM?=$(BLAST_DIR)/ublast.result.sam

USEARCH_INDEX?=$(BLAST_DIR)/usearch.index
USEARCH_INDEX_DONE?=$(BLAST_DIR)/.done_usearch_index
$(USEARCH_INDEX_DONE):
	$(call _start,$(BLAST_DIR))
	$(USEARCH_BIN) -makeudb_usearch $(BLAST_TARGET_FASTA) \
		-alpha $(USEARCH_TYPE) \
		-output $(USEARCH_INDEX)
	$(_end_touch)
ublast_index: $(USEARCH_INDEX_DONE)

USEARCH_DONE?=$(BLAST_DIR)/.done_usearch
$(USEARCH_DONE): $(USEARCH_INDEX_DONE)
	$(_start)
	$(call _time,$(BLAST_DIR),usearch) $(USEARCH_BIN) -usearch_global $(BLAST_QUERY_FASTA) \
		-db $(USEARCH_INDEX) \
		-threads $(USEARCH_THREADS) \
		-evalue $(USEARCH_EVALUE) \
		-id $(USEARCH_ID) \
		-strand both \
		-userout $(USEARCH_RESULT_SAM) -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits
	$(_end_touch)

USEARCH_PARSE_DONE?=$(BLAST_DIR)/.done_usearch_parse
$(USEARCH_PARSE_DONE): $(USEARCH_DONE)
	$(_start)
	perl $(_md)/pl/sam_parse.pl \
		$(USEARCH_RESULT_SAM) \
		$(BLAST_QUERY_TABLE) \
		$(BLAST_TARGET_TABLE) \
		nt \
		$(BLAST_RESULT)
	$(_end_touch)
ublast: $(USEARCH_PARSE_DONE)
