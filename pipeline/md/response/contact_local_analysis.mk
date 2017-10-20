LC_INPUT_DONE?=$(EE_BASE_MAP_DIR)/.done_input
$(LC_INPUT_DONE): $(RESPONSE_SELECT_NETWORK_DONE)
	$(_start)
	mkdir -p $(LC_DIR)
	perl $(_md)/pl/collect_associating_reads.pl \
		$(RESPONSE_NETWORK_SELECT) \
		$(RESPONSE_CONTIG_CLUSTER) \
		$(ANCHOR_TABLE) \
		$(PAIRED_DIR) \
		$(LC_READS)
	$(_end_touch)
lc_input: $(LC_INPUT_DONE)
