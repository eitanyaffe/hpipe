PP_COUNT_INPUT?=$(PREPROC_DIR)/.count_input
PP_COUNT_DUP?=$(PREPROC_DIR)/.count_dup
PP_COUNT_TRIM?=$(PREPROC_DIR)/.count_trim
PP_COUNT_ADAPT?=$(PREPROC_DIR)/.count_adaptor
PP_COUNT_LIGATION?=$(PREPROC_DIR)/.count_ligation_$(PREPROC_MODE)
PP_COUNT_HUMAN?=$(PREPROC_DIR)/.count_human_$(PREPROC_MODE)

##########################################################################################
# unite all source libraries and remove dups
##########################################################################################

FREE_PREFIX='*'
FREE_SUFFIX='*'

REMOVE_DUP_BIN?=$(_md)/bin.$(shell hostname)/remove_duplicates

INPUT_DIR=$(PREPROC_INTER_DIR)/input
INPUT_DIR_DONE=$(PREPROC_DIR)/.done_input
$(INPUT_DIR_DONE):
	$(call _start,$(INPUT_DIR))
	mkdir -p $(PREPROC_DIR)
	@echo "================================================================================"
	@echo Uniting all input fasta files and removing duplicates
	@echo "================================================================================"
	$(call _time,$(INPUT_DIR),input) perl $(_md)/pl/remove_duplicate_wrapper.pl \
		$(REMOVE_DUP_BIN) \
		$(INPUT_DIR)/R1 \
		$(INPUT_DIR)/R2 \
		$(LIB_COMPLEXITY_TABLE) \
		$(PP_COUNT_INPUT) \
		$(LIB_INPUT_DIRS)
	$(_end_touch)

$(PP_COUNT_DUP): $(INPUT_DIR_DONE)
	$(_md)/pl/count_fastq.pl $(INPUT_DIR) $(FREE_PREFIX) $(FREE_SUFFIX) duplicate $@

##########################################################################################
# split files
##########################################################################################

PP_SPLIT_DIR=$(PREPROC_INTER_DIR)/split

# split input fastq files
SPLIT_DONE=$(PREPROC_DIR)/.done_split
$(SPLIT_DONE): $(PP_COUNT_DUP)
	$(call _start,$(PP_SPLIT_DIR))
	$(_md)/pl/split_fastq.pl \
		$(PP_SPLIT_DIR) \
		$(SPLIT_SIZE) \
		F 0 0 \
		$(INPUT_DIR)/R1 \
		$(INPUT_DIR)/R2
ifeq ($(REMOVE_TRANSIENT),T)
	rm -rf $(INPUT_DIR)
endif
	$(_end_touch)

split: $(SPLIT_DONE)

##########################################################################################
# trim low quality (sickle)
##########################################################################################

SICKLE_QSUB_DIR?=$(QSUB_DIR)/preproc/$(LIB_ID)/sickle

TRIMMED_DIR=$(PREPROC_INTER_DIR)/trimmed
TRIMMED_DIR_DONE=$(PREPROC_DIR)/.done_trimmed
$(TRIMMED_DIR_DONE): $(SPLIT_DONE)
	$(call _start,$(TRIMMED_DIR))
	@echo "================================================================================"
	@echo "Trimming low quality ends (sickle)"
	@echo "================================================================================"
	$(call _time,$(PREPROC_DIR),trim) \
		$(_R) $(_md)/R/distrib_sickle.r distrib.sickle \
		idir=$(PP_SPLIT_DIR) \
		odir=$(TRIMMED_DIR) \
		sickle=$(SICKLE) \
		type=$(SICKLE_TYPE) \
		qsub.dir=$(SICKLE_QSUB_DIR) \
		batch.max.jobs=$(SICKLE_MAX_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		jobname=sickle
	@echo "================================================================================"
	$(_end_touch)
ifeq ($(REMOVE_TRANSIENT),T)
	rm -rf $(PP_SPLIT_DIR)
endif
trim_quality: $(TRIMMED_DIR_DONE)

$(PP_COUNT_TRIM): $(TRIMMED_DIR_DONE)
	$(_md)/pl/count_fastq.pl $(TRIMMED_DIR) $(FREE_PREFIX) $(FREE_SUFFIX) quality_trim $@

##########################################################################################
# remove adaptors (seqprep)
##########################################################################################

ADAPTOR_QSUB_DIR?=$(QSUB_DIR)/preproc/$(LIB_ID)/seqprep

NO_ADAPTOR_DIR=$(PREPROC_DIR)/no_adaptor
ADAPTOR_DONE?=$(PREPROC_DIR)/.done_adaptor
$(ADAPTOR_DONE): $(PP_COUNT_TRIM)
	$(call _start,$(NO_ADAPTOR_DIR))
	@echo "================================================================================"
	@echo "Removing adaptors (SeqPrep)"
	@echo "================================================================================"
	$(call _time,$(PREPROC_DIR),no_adaptor) \
		$(_R) $(_md)/R/distrib_remove_adaptors.r distrib.adpators \
		idir=$(TRIMMED_DIR) \
		odir=$(NO_ADAPTOR_DIR) \
		adaptor1=$(ADAPTOR1) adaptor2=$(ADAPTOR2) \
		seq.length=$(SEQPREP_LENGTH) \
		qsub.dir=$(ADAPTOR_QSUB_DIR) \
		batch.max.jobs=$(SEQPREP_MAX_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		seqprep=$(SEQPREP) \
		jobname=seqprep
	gunzip -f $(NO_ADAPTOR_DIR)/*.gz
	$(_end_touch)
ifeq ($(REMOVE_TRANSIENT),T)
	rm -rf $(TRIMMED_DIR) $(PREPROC_INTER_DIR)
endif

$(PP_COUNT_ADAPT): $(ADAPTOR_DONE)
	$(_md)/pl/count_fastq.pl $(NO_ADAPTOR_DIR) $(FREE_PREFIX) $(FREE_SUFFIX) remove_adaptors $@

##########################################################################################
# Remove Hi-C ligations (optional)
# Purpose: Use Hi-C library for genome assembly
##########################################################################################

SFILES=$(addprefix $(_md)/cpp/,clean_ligations.cpp Kmer.cpp Kmer.h util.cpp util.h)
$(eval $(call bin_rule2,clean_ligations,$(SFILES),-DKWORDS=$(CL_KWORDS)))
CLEAN_LIGATIONS_BIN=$(_md)/bin.$(shell hostname)/clean_ligations

CLEAN_LIGATIONS_DIR?=$(PREPROC_DIR)/clean_ligations_$(PREPROC_MODE)
CLEAN_LIGATIONS_DONE?=$(PREPROC_DIR)/.done_clean_$(PREPROC_MODE)
$(CLEAN_LIGATIONS_DONE): $(PP_COUNT_ADAPT)
	$(call _start,$(CLEAN_LIGATIONS_DIR))
ifeq ($(PREPROC_MODE),clean)
	$(call _time,$(PREPROC_DIR),no_ligations) $(CLEAN_LIGATIONS_BIN) \
		-mode split \
		-input_dir $(NO_ADAPTOR_DIR) \
		-input_side1 R1 -input_side2 R2 \
		-output_dir $(CLEAN_LIGATIONS_DIR) \
		-site $(CL_SITE) \
		-ksize $(CL_KSIZE) \
		-min_coverage $(CL_MIN_COVERAGE) \
		-min_ratio $(CL_MIN_RATIO) \
		-min_read_length $(CL_MIN_LENGTH) \
		-ofn_stats $(PREPROC_DIR)/clean_ligations_stats
endif
	$(_end_touch)

$(PP_COUNT_LIGATION): $(CLEAN_LIGATIONS_DONE)
ifeq ($(PREPROC_MODE),clean)
	$(_md)/pl/count_fastq.pl $(CLEAN_LIGATIONS_DIR) $(FREE_PREFIX) $(FREE_SUFFIX) remove_ligations $@
else
	touch $@
endif

##########################################################################################
# remove human dna (deconseq/bwa)
##########################################################################################

NO_HUMAN_DIR=$(PREPROC_DIR)/result_$(PREPROC_MODE)
NO_HUMAN_DONE=$(PREPROC_DIR)/.done_human_$(PREPROC_MODE)
NO_HUMAN_LOG=$(PREPROC_DIR)/.log_human_$(PREPROC_MODE)
NO_HUMAN_QSUB_DIR?=$(QSUB_DIR)/preproc/$(LIB_ID)/no_human

ifeq ($(PREPROC_MODE),clean)
NO_HUMAN_IDIR=$(CLEAN_LIGATIONS_DIR)
else
NO_HUMAN_IDIR=$(NO_ADAPTOR_DIR)
endif

HUMAN_SUFFIX='*clean*'

$(NO_HUMAN_DONE): $(PP_COUNT_LIGATION)
	$(call _start,$(NO_HUMAN_DIR))
ifeq ($(REMOVE_HUMAN),T)
	@echo "================================================================================"
	@echo "Removing human sequences (DeconSeq)"
	@echo "================================================================================"
	$(call _time,$(PREPROC_DIR),no_human) \
		$(_R) $(_md)/R/distrib_remove_human.r distrib.remove.human \
		deconseq=$(DECONSEQ_SCRIPT) \
		dbs=$(DECONSEQ_DBS) \
		wdir=$(DECONSEQ_DIR) \
		identity=$(DECONSEQ_IDENTITY) \
		coverage=$(DECONSEQ_COVERAGE) \
		idir=$(NO_HUMAN_IDIR) \
		odir=$(NO_HUMAN_DIR) \
		qsub.dir=$(NO_HUMAN_QSUB_DIR) \
		batch.max.jobs=$(DECONSEQ_MAX_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		jobname=deconseq
else
	$(foreach f,$(notdir $(wildcard $(NO_HUMAN_IDIR)/*.fastq)),cp $(NO_HUMAN_IDIR)/$f $(NO_HUMAN_DIR)/$(subst .fastq,,$f).clean.fq;$(ASSERT);)
endif
	$(_end_touch)
$(PP_COUNT_HUMAN): $(NO_HUMAN_DONE)
	$(_md)/pl/count_fastq.pl $(NO_HUMAN_DIR) $(FREE_PREFIX) $(HUMAN_SUFFIX) no_human_$(PREPROC_MODE) $@

##########################################################################################
# final links
##########################################################################################

PREPROC_FINAL_DONE?=$(PREPROC_FINAL_DIR)/.done
$(PREPROC_FINAL_DONE): $(PP_COUNT_HUMAN)
	$(call _start,$(PREPROC_FINAL_DIR))
	$(foreach f,$(notdir $(wildcard $(NO_HUMAN_DIR)/*clean.fq)),ln -sf ../../libs/$(LIB_ID)/result_$(PREPROC_MODE)/$f $(PREPROC_FINAL_DIR)/$f;$(ASSERT);)
	$(_end_touch)
preproc_result: $(PREPROC_FINAL_DONE)

##########################################################################################
# interface rule
##########################################################################################

make_preproc:
ifeq ($(PREPROC_MODES),simple)
	$(MAKE) preproc_result PREPROC_MODE=simple
else ifeq ($(PREPROC_MODES),clean)
	$(MAKE) preproc_result PREPROC_MODE=clean
else
	$(MAKE) preproc_result PREPROC_MODE=simple
	$(MAKE) preproc_result PREPROC_MODE=clean
endif

preproc_init: $(REMOVE_DUP_BIN) $(CLEAN_LIGATIONS_BIN)

.PHONY: make_preproc stats trim_adaptors split trim_quality input

force_remove_transient:
	rm -rf $(PREPROC_INTER_DIR)

########################################################################
# plots
########################################################################

plot_complex:
	$(_R) $(_md)/R/plots.r plot.complex \
		ldir=$(OUTDIR)/libs \
		fdir=$(SET_FIGURE_DIR)/lib_complexity \
		ids=$(PP_LIB_IDS) \
		titles=$(PP_LIB_TITLES) \
		cols=$(PP_COLS)

plot_pp_stats:
	$(_R) $(_md)/R/plots.r plot.stats \
		ldir=$(OUTDIR)/libs \
		fdir=$(SET_FIGURE_DIR)/lib_stats \
		ids=$(PP_LIB_IDS) \
		titles=$(PP_LIB_TITLES) \
		cols=$(PP_COLS)

make_plot_preproc: plot_complex plot_pp_stats
