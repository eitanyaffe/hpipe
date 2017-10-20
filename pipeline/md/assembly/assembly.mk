
RARIFY_DIR?=$(ASSEMBLY_DIR)/rarify
RARIFY_READS?=$(RARIFY_DIR)/reads.fasta

RARIFY_DONE?=$(ASSEMBLY_DIR)/.done_rarify
$(RARIFY_DONE):
	$(call _start,$(ASSEMBLY_DIR))
ifneq ($(RARIFY_FOLD),0)
	@mkdir -p $(RARIFY_DIR)
# broken here: script needs to be updated to support directories instead of pattern
	$(_md)/pl/rarify_reads.pl \
		$(RARIFY_FOLD) \
		$(RARIFY_READS) \
		$(ASSEMBLY_INPUT_NAME_PATTERN) \
		$(ASSEMBLY_INPUT_DIRS)
endif
	$(_end_touch)

ASSEMBLY_INPUT_FILE_TABLE?=$(ASSEMBLY_DIR)/input_table
$(ASSEMBLY_INPUT_FILE_TABLE): $(RARIFY_DONE)
	$(call _start,$(ASSEMBLY_DIR))
ifeq ($(RARIFY_FOLD),0)
	$(call _assert,ASSEMBLY_LIB_IDS)
	find $(ASSEMBLY_INPUT_DIRS) -name "$(ASSEMBLY_INPUT_NAME_PATTERN)" > $@
else
	find $(RARIFY_READS) > $@
endif

include $(_md)/assembly_$(ASSEMBLER).mk

CONTIG_GC_TABLE_DONE=$(ASSEMBLY_DIR)/.done_gc
$(CONTIG_GC_TABLE_DONE): $(FULL_CONTIG_FILE)
	$(_start)
	$(_md)/pl/compute_gc.pl $(FULL_CONTIG_FILE) $(CONTIG_GC_BINSIZE) $(CONTIG_GC_BINNED) $(CONTIG_GC_TABLE)
	$(_end_touch)
assembly_gc: $(CONTIG_GC_TABLE_DONE)

$(CONTIG_FILE): $(FULL_CONTIG_TABLE)
	$(_start)
	$(_R) $(_md)/R/select_contigs.r top \
	        table=$(FULL_CONTIG_TABLE) \
		ofn=$(CONTIG_TABLE) \
		min.length=$(ASSEMBLY_MIN_LEN)
	$(_md)/pl/select_contigs.pl \
		$(FULL_CONTIG_FILE) \
		$(CONTIG_TABLE) \
		$@
	$(_end)
make_assembly: $(CONTIG_FILE) $(CONTIG_GC_TABLE_DONE)
	@echo "assembly done, output is in directory ASSEMBLY_DIR=$(ASSEMBLY_DIR)"

make_plot_assembly:
	$(_R) R/assembly.r plot.assembly \
		ids=$(ASSEMBLY_IDS) \
		titles=$(ASSEMBLY_TITLES) \
		cols=$(ASSEMBLY_COLS) \
		idir=$(OUTDIR)/assembly \
		min.length=$(ASSEMBLY_MIN_LEN) \
		fdir=$(SET_FIGURE_DIR)/assembly

###############################################################################
# make rarified assemblies
###############################################################################

rarify_assembly:
	$(foreach v,$(RARIFY_FOLDS),$(MAKE) make_assembly RARIFY_FOLD=$(v); $(ASSERT);)

$(RARIFY_TABLE):
	@$(MAKE) rarify_assembly
	$(_R) $(_md)/R/rarify_reads.r compute.rarify \
		folds=$(RARIFY_FOLDS) \
		base.dir=$(ASSEMBLY_BASEDIR) \
		fold0.table=$(ASSEMBLY_INPUT_FILE_TABLE) \
		ofn=$@
rarify: $(RARIFY_TABLE)

plot_rarify: $(RARIFY_TABLE)
	$(_R) $(_md)/R/rarify_reads.r plot.rarify \
		ifn=$(RARIFY_TABLE) \
		fdir=$(FIGURE_DIR)/rarify
