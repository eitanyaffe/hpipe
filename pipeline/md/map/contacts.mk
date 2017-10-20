# cross coverage
CROSS_COVERAGE_DONE?=$(CROSS_COVERAGE_DIR)/.done
$(CROSS_COVERAGE_DONE): $(PAIRED_DONE)
	$(call _start,$(CROSS_COVERAGE_DIR))
	$(call TIME,$(CROSS_COVERAGE_DIR)) $(_md)/pl/cross_coverage.pl \
		$(CONTIG_TABLE) \
		$(PAIRED_DIR) \
		$(CROSS_BINSIZE) \
		$(CROSS_MIN_DIST) \
		$(CROSS_MAX_DIST) \
		$(CROSS_COVERAGE_DIR)
	$(_end_touch)
cross: $(CROSS_COVERAGE_DONE)

INTRA_CONTIG_CONTACT_DONE?=$(INTRA_CONTIG_CONTACT_DIR)/.done
$(INTRA_CONTIG_CONTACT_DONE):
	$(call _start,$(INTRA_CONTIG_CONTACT_DIR))
	$(call TIME,$(INTRA_CONTIG_CONTACT_DIR)) $(_md)/pl/intra_contig_contacts.pl \
		$(CONTIG_TABLE) \
		$(PAIRED_DIR) \
		$(INTRA_CONTIG_CONTACT_DIR)
	$(_end_touch)
intra_contig: $(INTRA_CONTIG_CONTACT_DONE)

intra_contig_plot: $(INTRA_CONTIG_CONTACT_DONE)
	$(_start)
	$(_R) $(_md)/R/intra_contig.r plot.intra \
		ifn=$(CONTIG_TABLE) \
		idir.coverage=$(COVERAGE_DIR) \
		idir.contacts=$(INTRA_CONTIG_CONTACT_DIR) \
		threshold=$(INTRA_CONTIG_MIN_LENGTH) \
		odir=$(INTRA_FDIR)
	$(_end)

.PHONY: cross
