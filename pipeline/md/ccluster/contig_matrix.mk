CONTIG_MATRIX_DONE=$(CONTIG_MATRIX_DIR)/.done
$(CONTIG_MATRIX_DONE): $(SIM_DONE)
	$(call _start,$(CONTIG_MATRIX_DIR))
	@mkdir -p \
		$(CONTIG_CONTACTS) \
		$(CONTIG_MASKED_CONTACTS)
	$(call _time,$(DATASET_DIR),matrix) \
		perl $(_md)/pl/contig_matrix.pl \
			$(PAIRED_DIR) \
			$(CONTIG_TABLE) \
			$(COVERAGE_TABLE) \
			$(SIM_RESULT_DIR) \
			$(CONTIG_MATRIX_SIMILARITY_OFFSET) \
			$(CONTIG_CONTACTS) \
			$(CONTIG_MASKED_CONTACTS) \
			$(CONTIG_MATRIX)
	$(_end_touch)

CONTIG_MATRIX_FILTERED_DONE?=$(CONTIG_MATRIX_DIR)/.done_filter
$(CONTIG_MATRIX_FILTERED_DONE): $(CONTIG_MATRIX_DONE)
	$(_start)
	$(_R) R/contig_matrix.r filter.noise \
		ifn=$(CONTIG_MATRIX) \
		ofn.filter.params=$(CONTIG_MATRIX_FILTER_PARAMS) \
		ofn.stats=$(CONTIG_MATRIX_STATS) \
		ofn=$(CONTIG_MATRIX_FILTERED)
	$(_end_touch)
filter_noise: $(CONTIG_MATRIX_FILTERED_DONE)

contig_matrix: $(CONTIG_MATRIX_FILTERED_DONE)

make_plot_contig_matrix:
	$(_R) R/contig_matrix.r plot.contig.matrix \
		adir=$(ASSEMBLY_DIR) \
		fdir=$(SET_FIGURE_DIR)/cc_matrix \
		aid=$(ASSEMBLY_ID) \
		ids=$(LIB_IDS) \
		titles=$(LIB_TITLES)

plot_filter_noise:
	$(_R) R/contig_matrix.r plot.filter.noise \
		ifn=$(CONTIG_MATRIX) \
		ifn.mdl=$(CONTIG_MATRIX_FILTER_PARAMS) \
		fdir=$(CCLUSTER_FIGURE_DIR)

.PHONY: contig_matrix plot_filter_noise filter_noise
