########################################################
# Input
# MDL_FENDS: input fends
# MDL_SCOPE: [inter_anchor|intra_anchor] infer model over limited fend pairs
# MDL_MAT: input fend pairs
# MDL_ANCHORS: table of contig clustering into anchors
# MDL_ANCHOR_KEY/VALUE: fields in anchor table
# MDL_MFN: file with model fields
#
# Output
# MDL_DIR: directory
# MDL_QSUB: directory to use for temp files
########################################################

MDL_PREFIX?=$(MDL_DIR)/main
MDL_FENDS_ANCHORS?=$(MDL_PREFIX).anchors
MDL_FENDS_ABUNDANCE?=$(MDL_PREFIX).abundance
MDL_FENDS_BINNED?=$(MDL_PREFIX).binned
MDL_FENDS_RESTRICTED?=$(MDL_PREFIX).restricted
MDL_ABUN_INIT_FUNC?=$(MDL_PREFIX)_abundance_bin.f

MDL_COUNTS?=$(MDL_PREFIX).counts
MDL_TOTAL_COUNTS?=$(MDL_PREFIX).total_counts
MDL_NM?=$(MDL_PREFIX).nm

MDL_TMP_FILES_DIR?=$(MDL_QSUB)/files
MDL_TMP_SCRIPTS_DIR?=$(MDL_QSUB)/scripts

#############################################################################################
# generating fend table with anchors and model bins
#############################################################################################

$(MDL_FENDS_ANCHORS):
	$(call _start,$(MDL_DIR))
	$(_md)/pl/assign_anchor.pl \
		$(MDL_ANCHORS) \
		$(MDL_ANCHOR_KEY) \
		$(MDL_ANCHOR_VALUE) \
		$(MDL_FENDS) \
		$@
	$(_end)

$(MDL_FENDS_ABUNDANCE): $(MDL_FENDS_ANCHORS)
	$(_start)
	$(_R) R/bin_abundance.r bin.abundance \
	           fends.ifn=$(MDL_FENDS_ANCHORS) \
	           mfn=$(MDL_MFN) \
	           bin.ranges.ofn=$(MDL_PREFIX)_abundance.bin_ranges \
	           correction.ofn=$(MDL_ABUN_INIT_FUNC) \
	           fends.ofn=$@
	$(_end)
model_bin_abundance: $(MDL_FENDS_ABUNDANCE)

$(MDL_FENDS_BINNED): $(MDL_FENDS_ABUNDANCE)
	$(_start)
	$(_R) R/bin_fields.r bin.fields.mfn \
		fends.fn=$(MDL_FENDS_ABUNDANCE) \
		ofn.prefix=$(MDL_PREFIX) \
		mfn=$(MDL_MFN) \
		round=3
	$(_end)

DONE_BIN_FENDS?=$(MDL_DIR)/.done_binned_fends
$(DONE_BIN_FENDS): $(MDL_FENDS_BINNED)
	$(_start)
	$(_md)/pl/anchor_only_fends.pl \
		$(MDL_FENDS_BINNED) \
		anchor \
		$(MDL_FENDS_RESTRICTED)
	$(_end_touch)
model_fends: $(DONE_BIN_FENDS)

#############################################################################################
# observed and total counts (nm)
#############################################################################################

$(MDL_COUNTS): $(DONE_BIN_FENDS)
	$(_start)
	cp $(MDL_MAT) $(MDL_PREFIX).mat
	$(_md)/pl/compute_n_mfn.pl \
		$(MDL_PREFIX) \
		$(MDL_SCOPE) \
		F restricted \
		$(MDL_MFN) \
		$(_md)
	$(_end)
model_counts: $(MDL_COUNTS)

$(MDL_TOTAL_COUNTS): $(MDL_FENDS_BINNED)
	$(_start)
	$(_R) R/model_predict.r compute.total.counts \
		prefix=$(MDL_PREFIX) \
		model.ifn=$(MDL_MFN) \
		tmp.dir=$(MDL_TMP_FILES_DIR) \
		log.dir=$(MDL_TMP_SCRIPTS_DIR) \
		binary=$(MDL_BINARY) \
		wd=$(_md) \
		dtype=$(DTYPE) \
		max.jobs.fn=$(MAX_JOBS_FN) \
		fends.suffix=restricted \
		scope=$(MDL_SCOPE) \
		Rcall="$(_Rcall)"
	$(_end)
model_total_counts: $(MDL_TOTAL_COUNTS)

$(MDL_NM): $(MDL_COUNTS) $(MDL_TOTAL_COUNTS)
	$(_start)
	$(_md)/pl/append_table_mfn.pl \
		$(MDL_TOTAL_COUNTS) \
		$(MDL_COUNTS) \
                $@ \
		count 0 F \
		$(MDL_MFN) \
		$(_md)
	$(_end)
model_nm: $(MDL_NM)

#############################################################################################
# model
#############################################################################################

MDL_DONE?=$(MDL_DIR)/.done
$(MDL_DONE): $(MDL_NM)
	$(_start)
	$(_R) R/model.r learn.model \
		model.prefix=$(MDL_PREFIX) \
		model.fn=$(MDL_MFN)
	$(_end_touch)
model: $(MDL_DONE)

#############################################################################################
# gcc compilation
#############################################################################################

model_init_bin: $(MDL_BINARY)

#############################################################################################
# plotting
#############################################################################################

# bias model matrices
#BIAS_BREAKS?="-8 -4 0 4 8"
#BIAS_BREAKS?="-7 -6.2 -5.5 -4.8 -4"
#BIAS_BREAKS?="-7 0 5 -10 20"
#BIAS_BREAKS?="0 5 10 15 20"
BIAS_COLORS?="black blue white orange red"
BIAS_BREAKS?="-5 0 5 7.5 10"
plot_model:
	$(_start)
	$(_R) R/plot_model_matrices.r plot.feature.matrices \
	  	idir=$(MDL_DIR) \
		dataset=$(DATASET) \
		fdir=$(MDL_FDIR) \
	  	model.fn=$(MDL_MFN) \
		wd=$(CURDIR) \
		width=6 height=7 \
	  	breaks=$(BIAS_BREAKS) \
		break.cols=$(BIAS_COLORS)
	$(_end)

