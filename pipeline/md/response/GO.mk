####################################################################
# GO over network genes
####################################################################

RESPONSE_GO_CTRL_DONE?=$(EE_BASE_MAP_DIR)/.done_GO_ctrl
$(RESPONSE_GO_CTRL_DONE): $(RESPONSE_SELECT_NETWORK_DONE)
	$(_start)
	$(_md)/pl/GO_append.pl \
		$(UNIREF_GENE_TAX_TABLE) \
		$(UNIREF_GENE_GO) \
		$(GO_TREE) \
		$(RESPONSE_GO_TABLE_CTRL)
	$(_md)/pl/GO_analysis.pl \
		$(RESPONSE_GO_TABLE_CTRL) \
		$(RESPONSE_GO_MIN_IDENTITY) \
		$(GO_TREE) \
		$(RESPONSE_GO_SUMMARY_CTRL_PREFIX)
	$(_end_touch)
element_GO_ctrl: $(RESPONSE_GO_CTRL_DONE)

RESPONSE_GO_DONE?=$(EE_BASE_MAP_DIR)/.done_GO
$(RESPONSE_GO_DONE): $(RESPONSE_SELECT_NETWORK_DONE)
	$(_start)
	$(_md)/pl/GO_append.pl \
		$(RESPONSE_NETWORK_GENES_SELECT) \
		$(UNIREF_GENE_GO) \
		$(GO_TREE) \
		$(RESPONSE_GO_TABLE)
	$(_md)/pl/GO_analysis.pl \
		$(RESPONSE_GO_TABLE) \
		$(RESPONSE_GO_MIN_IDENTITY) \
		$(GO_TREE) \
		$(RESPONSE_GO_SUMMARY_PREFIX)
	$(_end_touch)
element_GO: $(RESPONSE_GO_DONE)

MERGE_RESPONSE_GO_DONE?=$(EE_BASE_MAP_DIR)/.done_merge_GO
$(MERGE_RESPONSE_GO_DONE): $(RESPONSE_GO_DONE) $(RESPONSE_GO_CTRL_DONE)
	$(_start)
	$(_R) R/merge_GO.r merge.GO \
		ifn.genes=$(RESPONSE_NETWORK_GENES_SELECT) \
		ifn.genes.ctrl=$(UNIREF_GENE_TAX_TABLE) \
		ifn.prefix=$(RESPONSE_GO_SUMMARY_PREFIX) \
		ifn.prefix.ctrl=$(RESPONSE_GO_SUMMARY_CTRL_PREFIX) \
		ofn=$(RESPONSE_GO_MERGE)
	$(_end_touch)
merge_GO: $(MERGE_RESPONSE_GO_DONE)
