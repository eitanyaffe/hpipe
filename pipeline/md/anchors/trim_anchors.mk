##############################################################################
# AMAP_IN_FENDS: fends table
# AMAP_IN_MAT: fend pairs
# MODEL_ID: string id of model
# AMAP_IN_MODEL_PREFIX: filename prefix of model files
# AMAP_OUT_DIR: output directory
##############################################################################

AMAP_OUT_BINNED?=$(AMAP_OUT_DIR)/fends.binned
AMAP_OUT_CONTIG_BINS?=$(AMAP_OUT_DIR)/contig.bins
AMAP_OUT_OBSERVED?=$(AMAP_OUT_DIR)/observed
AMAP_OUT_EXPECTED?=$(AMAP_OUT_DIR)/expected
AMAP_OUT_RESULT?=$(AMAP_OUT_DIR)/united

AMAP_POST_TRIM_MATRIX?=$(AMAP_OUT_DIR)/multi_clean.matrix
AMAP_TRIM_TABLE?=$(AMAP_OUT_DIR)/contigs.trimmed

AMAP_TMP_DIR?=$(QSUB_ANCHOR_DIR)/trim_anchor/scripts
AMAP_LOG_DIR?=$(QSUB_ANCHOR_DIR)/trim_anchor/log

#AMAP_OUT_FDIR?=$(MAP_FIGURE_DIR)/anchors/$(ANCHOR)

AMAP_OUT_BINNED_DONE?=$(AMAP_OUT_DIR)/.done_binning_contigs
$(AMAP_OUT_BINNED_DONE):
	$(call _start,$(AMAP_OUT_DIR))
	$(_md)/pl/bin_contigs.pl \
		$(AMAP_IN_FENDS) contig \
		$(AMAP_OUT_BINNED) \
		$(AMAP_OUT_CONTIG_BINS)
	$(_end_touch)

$(AMAP_OUT_OBSERVED): $(AMAP_OUT_BINNED_DONE)
	$(_start)
	$(call _time,$(AMAP_OUT_DIR),observed) \
		$(_md)/pl/cc_observed.pl \
			$(AMAP_OUT_BINNED) anchor \
			$(AMAP_IN_MAT) \
			$(AMAP_MULTI_BIN_SIZE) \
			$@
	$(end)
trim_obs: $(AMAP_OUT_OBSERVED)

$(AMAP_OUT_EXPECTED): $(AMAP_OUT_BINNED_DONE)
	$(_start)
	@rm -rf $(AMAP_TMP_DIR) $(AMAP_LOG_DIR)
	$(call _time,$(AMAP_OUT_DIR),expected) \
		$(_R) R/model_predict.r compute.expected.counts \
		   binary=$(MODEL_INTEGRATE_BINARY) \
		   tmp.dir=$(AMAP_TMP_DIR) \
		   log.dir=$(AMAP_LOG_DIR) \
		   wd=$(_md) \
	           model.prefix=$(AMAP_IN_MODEL_PREFIX) \
	           fends.ifn=$(AMAP_OUT_BINNED) \
	  	   model.ifn=$(AMAP_MFN) \
		   scope=anchored model=std \
	           ofields.x=contig_bin ofields.y=anchor omit.y.zero=T \
	  	   num.splits=$(AMAP_NUM_SPLIT) \
		   max.jobs.fn=$(MAX_JOBS_FN) \
		   req.mem=$(AMAP_REQ_MEM) \
		   dtype=$(DTYPE) \
		   Rcall="$(_Rcall)" \
	    	   ofn=$@
	$(_end)
trim_exp: $(AMAP_OUT_EXPECTED)

#############################################################################################################################
# compute cell-contig obs and expected matrix
#############################################################################################################################

$(AMAP_OUT_RESULT): $(AMAP_OUT_OBSERVED) $(AMAP_OUT_EXPECTED)
	$(_start)
	$(_md)/pl/trim_unite.pl \
		$(AMAP_OUT_CONTIG_BINS) \
		$(INITIAL_ANCHOR_TABLE) cluster \
		$(AMAP_OUT_OBSERVED) T \
		$(AMAP_OUT_EXPECTED) \
		$@
	$(_end)
trim_unite: $(AMAP_OUT_RESULT)

# only contigs cells associated with cells
AMAP_OUT_FILTERED?=$(AMAP_OUT_DIR)/united.filtered

$(AMAP_OUT_FILTERED): $(AMAP_OUT_RESULT)
	$(_start)
	$(_md)/pl/filter_no_cell_contigs.pl \
		$(AMAP_OUT_RESULT) \
		$@
	$(_end)

trim_map: $(AMAP_OUT_FILTERED)

#############################################################################################################################
# remove multi contigs, optimized code
#############################################################################################################################

AMAP_MULTI_THRESHOLD_TABLE?=$(AMAP_OUT_DIR)/threshold_table
$(AMAP_MULTI_THRESHOLD_TABLE): $(AMAP_OUT_FILTERED)
	$(_start)
	$(_R) R/trim_anchors.r compute.trim.threshold \
		ifn=$(AMAP_OUT_FILTERED) \
		min.contacts=$(AMAP_MULTI_MIN_CONTACTS) \
		min.enrichment=$(AMAP_MULTI_MIN_ENRICHMENT) \
		fdr=$(AMAP_TRIM_FDR) \
		ofn=$@
	$(_end)
ttable: $(AMAP_MULTI_THRESHOLD_TABLE)

plot_trim_in: $(AMAP_MULTI_THRESHOLD_TABLE)
	$(_start)
	$(_R) R/trim_anchors.r plot.trim.input \
		ifn=$(AMAP_OUT_FILTERED) \
		min.contacts=$(AMAP_MULTI_MIN_CONTACTS) \
		ifn.threshold=$(AMAP_MULTI_THRESHOLD_TABLE) \
		contigs.fn=$(CONTIG_TABLE) \
		fdir=$(ANCHOR_FIGURE_DIR)/trim_anchor/in
	$(_end)

# debug: compare model of model_integrate and model of trim_anchors
AMAP_VERIFY_CA_MATRIX?=F

FILES=$(addprefix $(_md)/cpp/,Model.cpp Model.h Params.cpp Params.h util.h util.cpp)
$(eval $(call bin_rule2,trim_anchors,$(FILES)))
TRIM_ANCHORS_BINARY=$(_md)/bin.$(shell hostname)/trim_anchors

MDL_PARAMS=$(shell Rscript $(_md)/R/model_parse.r $(INIT_MDL_MFN) $(AMAP_IN_MODEL_PREFIX))

TRIM_ANCHOR_DONE?=$(AMAP_OUT_DIR)/.done_trim
TRIM_ANCHOR_LOG_OUT?=$(AMAP_OUT_DIR)/.log_trim_stdout
TRIM_ANCHOR_LOG_ERR?=$(AMAP_OUT_DIR)/.log_trim_stderr
$(TRIM_ANCHOR_DONE): $(AMAP_MULTI_THRESHOLD_TABLE)
	$(_start)
	$(call _time,$(AMAP_OUT_DIR),trim) $(TRIM_ANCHORS_BINARY) \
		-fends $(AMAP_OUT_BINNED) \
		-contacts $(AMAP_IN_MAT) \
		-contigs $(CONTIG_TABLE) \
		-threshold_table $(AMAP_MULTI_THRESHOLD_TABLE) \
		-min_contacts $(AMAP_MULTI_MIN_CONTACTS) \
		-min_anchor_size $(AMAP_MIN_LENGTH) \
		-ca_matrix $(AMAP_OUT_RESULT) \
		-verify_ca_matrix $(AMAP_VERIFY_CA_MATRIX) \
		-prior_fn $(AMAP_IN_MODEL_PREFIX).prior \
		$(MDL_PARAMS) \
		-model_num 2 \
		-ofn_matrix $(AMAP_POST_TRIM_MATRIX) \
		-ofn $(AMAP_TRIM_TABLE) \
		 > $(TRIM_ANCHOR_LOG_OUT) 2> $(TRIM_ANCHOR_LOG_ERR)
	$(_end_touch)

#############################################################################################################################
# plot anchor coverage
#############################################################################################################################

plot_anchor_coverage:
	$(_R) R/trim_anchors.r plot.anchor.coverage \
		ifn=$(AMAP_TRIM_TABLE) \
		ifn.coverage=$(COVERAGE_TABLE) \
		factor=$(AMAP_TRIM_COVERAGE_FACTOR) \
		fdir=$(ANCHOR_FIGURE_DIR)/trim_anchor/coverage

#############################################################################################################################
# select final anchors, omitting coverage outliers
#############################################################################################################################

TRIM_ANCHOR_SELECT_DONE?=$(AMAP_OUT_DIR)/.done_trim_select
$(TRIM_ANCHOR_SELECT_DONE): $(TRIM_ANCHOR_DONE)
	$(call _start,$(ANCHOR_DIR))
	$(_R) R/trim_anchors.r select.contigs \
		ifn=$(AMAP_TRIM_TABLE) \
		ifn.coverage=$(COVERAGE_TABLE) \
		factor=$(AMAP_TRIM_COVERAGE_FACTOR) \
		max.sd=$(AMAP_TRIM_COVERAGE_MAX_SD) \
		ofn.contigs=$(ANCHOR_TABLE) \
		ofn.lookup=$(ANCHOR_RENAME_LOOKUP_TABLE)
	$(_end_touch)

plot_trim_out: $(TRIM_ANCHOR_SELECT_DONE)
	$(_start)
	$(_R) R/trim_anchors.r plot.trim \
		pre.cc=$(AMAP_OUT_RESULT) \
		post.cc=$(AMAP_POST_TRIM_MATRIX) \
		ifn.contigs=$(AMAP_TRIM_TABLE) \
		min.contacts=$(AMAP_MULTI_MIN_CONTACTS) \
		ifn.threshold=$(AMAP_MULTI_THRESHOLD_TABLE) \
		fdir=$(ANCHOR_FIGURE_DIR)/trim_anchor/out
	$(_end)

#############################################################################################################################
#############################################################################################################################


trim_anchors: $(TRIM_ANCHOR_SELECT_DONE)
	@echo "DONE Anchors, ANCHOR_TABLE=$(ANCHOR_TABLE)"

plot_trim: plot_anchor_coverage plot_trim_out

plot_trim_extended: plot_trim_in

trim_init: $(TRIM_ANCHORS_BINARY)

$(call _add_module_target,$m,trim_anchors,infer anchor contig sets)
$(call _add_module_target,$m,trim_plots,plot final anchor figures)

