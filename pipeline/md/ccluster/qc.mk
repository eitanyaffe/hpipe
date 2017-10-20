
# contig-based analysis
$(MAP_CONTACT_TABLE): $(CONTIG_MATRIX_DONE)
	$(_start)
	$(_md)/pl/contig_analysis.pl \
		$(CONTIG_TABLE) \
		$(CONTIG_CONTACTS) \
		$(MAP_CONTACT_CLOSE_CIS_THRESHOLD) \
		$(MAP_CONTACT_FAR_CIS_THRESHOLD) \
		$(MAP_CONTACT_BINSIZE) \
		$(MAP_CONTACT_TABLE_MAX_READS) \
		$@
	$(_end)
contig_analysis: $(MAP_CONTACT_TABLE)

plot_contig_analysis:
	$(_R) $(_md)/R/contact_analysis.r plot.info \
		adir=$(ASSEMBLY_DIR) \
		mreads=$(MAP_CONTACT_TABLE_MAX_READS) \
		fdir=$(SET_FIGURE_DIR)/contact_analysis/mreads_$(MAP_CONTACT_TABLE_MAX_READS) \
		aid=$(ASSEMBLY_ID) \
		ids=$(LIB_IDS) \
		titles=$(LIB_TITLES)

# cis-decay summary analysis
CIS_DECAY_DONE=$(CIS_DECAY_DIR)/.done
$(CIS_DECAY_DONE): $(CONTIG_MATRIX_DONE)
	$(call _start,$(CIS_DECAY_DIR))
	$(_md)/pl/cis_decay.pl \
		$(CONTIG_TABLE) \
		$(CONTIG_CONTACTS) \
		$(CIS_DECAY_BIN_LOG_START) \
		$(CIS_DECAY_BIN_LOG_END) \
		$(CIS_DECAY_BIN_LOG_STEP) \
		$(CIS_DECAY_GAP) \
		$(CIS_DECAY_MAX_READS) \
		$(CIS_DECAY_TABLE) \
		$(CIS_DECAY_SUMMARY)
	$(_end_touch)
make_cis_decay: $(CIS_DECAY_DONE)

make_plot_cis_decay:
	$(_R) $(_md)/R/cis_decay.r plot.cis.decay \
		ifn=$(CIS_DECAY_TABLE) \
		id=$(LIB_ID) \
		distal.threshold=$(CIS_DECAY_DISTAL_THRESHOLD) \
		fdir=$(SET_FIGURE_DIR)/cis_decay/$(SET_TITLE)/$(LIB_ID) \
		ylim.n=$(CIS_DECAY_YLIM_N)

# assembly summary
make_plot_cis_decay_assembly:
	$(_R) $(_md)/R/cis_decay.r plot.summary \
		distal.threshold=$(CIS_DECAY_DISTAL_THRESHOLD) \
		adir=$(ASSEMBLY_DIR) \
		mreads=$(MAP_CONTACT_TABLE_MAX_READS) \
		fdir=$(SET_FIGURE_DIR)/cis_decay/$(SET_TITLE)/summary \
		aid=$(ASSEMBLY_ID) \
		ids=$(LIB_IDS) \
		titles=$(LIB_TITLES) \
		ymax=$(CIS_DECAY_SUMMARY_YMAX)

# compare contig matrix between replicates
compare_cc_matrix:
	$(_R) $(_md)/R/plot_compare_cc_matrix.r plot.compare.cc.matrix \
		adir=$(ASSEMBLY_DIR) \
		fdir=$(SET_FIGURE_DIR)/compare_cc_matrix \
		aid=$(ASSEMBLY_ID) \
		ids=$(LIB_IDS) \
		titles=$(LIB_TITLES)
