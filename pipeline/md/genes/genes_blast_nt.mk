
BLAST_BIN?=blastn

# params
BLAST_THREADS=40
BLAST_EVALUE=0.001

# files
BLAST_RESULT_SAM?=$(BLAST_DIR)/blast.result.sam

BLAST_INDEX_NT_DONE?=$(BLAST_TARGET_FASTA).done_index
$(BLAST_INDEX_NT_DONE):
	$(call _start,$(BLAST_DIR))
	$(call _time,$(BLAST_DIR),index) makeblastdb \
		-in $(BLAST_TARGET_FASTA) \
		-parse_seqids -dbtype nucl
	$(_end_touch)
blast_index: $(BLAST_INDEX_NT_DONE)

BLAST_NT_DONE?=$(BLAST_DIR)/.done_raw_nt
$(BLAST_NT_DONE): $(BLAST_INDEX_NT_DONE)
	$(call _start,$(BLAST_DIR))
	$(call _time,$(BLAST_DIR),blast) blastn \
		-query $(BLAST_QUERY_FASTA) \
		-task blastn \
		-db $(BLAST_QUERY_FASTA) \
		-out $(BLAST_RESULT_SAM) \
		-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
		-evalue $(BLAST_EVALUE) \
		-num_threads $(BLAST_THREADS)
	$(_end_touch)

BLAST_PARSE_NT_DONE?=$(BLAST_DIR)/.done_parse_nt
$(BLAST_PARSE_NT_DONE): $(BLAST_NT_DONE)
	$(_start)
	$(call _time,$(BLAST_DIR),parse_sam) perl $(_md)/pl/sam_parse.pl \
		$(BLAST_RESULT_SAM) \
		$(BLAST_QUERY_TABLE) \
		$(BLAST_TARGET_TABLE) \
		nt \
		$(BLAST_RESULT)
	$(_end_touch)
blast_nt: $(BLAST_PARSE_NT_DONE)
