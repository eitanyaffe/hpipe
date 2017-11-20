
export_table:
	$(call _start,$(ANCHOR_EXPORT_DIR))
	$(_R) R/export.r export.table \
		ofn=$(EXPORT_TABLE) \
		$(call _export_variable,$(EXPORT_VARIALBES)) \
		$(EXPORT_VARIALBES_NOEVAL)
	$(_end)

export:
	$(call _start,$(EXPORT_ODIR))
	$(_R) R/export.r export \
		odir=$(EXPORT_ODIR) \
		$(call _export_variable,$(EXPORT_VARIALBES)) \
		$(EXPORT_VARIALBES_NOEVAL)
	tar cvf $(EXPORT_ODIR_TAR) $(EXPORT_ODIR)
	$(_end)

flat_basic:
	$(call _start,$(FLAT_DIR))
	cp $(FULL_CONTIG_TABLE) $(FLAT_DIR)/contig.table
	cp $(CONTIG_FILE) $(FLAT_DIR)/contig.fasta
	$(_R) R/export.r copy.table ifn=$(ANCHOR_TABLE) ofn=$(FLAT_DIR)/anchor.table \
		fields.in="contig anchor" \
		fields.out="contig anchor"
	$(_R) R/export.r copy.table ifn=$(CA_ANCHOR_CONTIGS) ofn=$(FLAT_DIR)/contig_anchor.table \
		fields.in="contig anchor contig_total_count contig_expected enrichment" \
		fields.out="contig anchor observed expected score"
	mkdir -p $(FLAT_DIR)/model
	cp $(FINAL_MODEL_PREFIX).binned $(FLAT_DIR)/model/fends.table
	cp $(FINAL_MODEL_PREFIX)_frag_len_bin.f $(FLAT_DIR)/model/fragment_length.f
	cp $(FINAL_MODEL_PREFIX)_frag_len.bin_ranges $(FLAT_DIR)/model/fragment_length.bins
	cp $(FINAL_MODEL_PREFIX)_abundance_bin.f $(FLAT_DIR)/model/abundance.f
	cp $(FINAL_MODEL_PREFIX)_abundance.bin_ranges $(FLAT_DIR)/model/abundance.bins
	cp $(FINAL_MODEL_PREFIX).prior $(FLAT_DIR)/model/prior
	$(_end)
