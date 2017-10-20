# fends
EE_MAP_FENDS?=$(EE_MAP_DIR)/fends
EE_MAP_BINNED?=$(EE_MAP_DIR)/fends.binned

# observed
EE_MAP_OBSERVED?=$(EE_MAP_DIR)/count.o

# expected
EE_MAP_EXPECTED?=$(EE_MAP_DIR)/count.e

#######################################################################################################################
# prepare fends table
########################################################################################################################

$(EE_MAP_FENDS):
	$(call _start,$(EE_MAP_DIR))
	$(_md)/pl/ee_fends.pl \
		$(EE_MAP_IN_FENDS) \
		$(EE_MAP_IN_TABLE) \
		$@
	$(_end)
assign: $(EE_MAP_FENDS)

$(EE_MAP_BINNED): $(EE_MAP_FENDS)
	$(_start)
	$(_R) R/bin_fields_permute.r bin.fields \
		ifn=$(EE_MAP_FENDS) \
		ofn=$(EE_MAP_BINNED) \
		permute=T \
		field=cluster
	$(_end)
ee_binned: $(EE_MAP_BINNED)

########################################################################################################################
# observed
########################################################################################################################

$(EE_MAP_OBSERVED): $(EE_MAP_BINNED)
	$(_start)
	$(_md)/pl/ee_observed.pl \
		$(EE_MAP_BINNED) \
		$(EE_MAP_IN_MAT) \
		$@
	$(_end)
ee_obs: $(EE_MAP_OBSERVED)

########################################################################################################################
# expected contigs
########################################################################################################################

EE_QSUB_DIR?=$(QSUB_DIR)/datasets/$(DATASET)/ee_map

EE_TMP_CONTIG_DIR?=$(EE_QSUB_DIR)/contig_tmp
EE_LOG_CONTIG_DIR?=$(EE_QSUB_DIR)/contig_log

# MODEL_INTEGRATE_BINARY?=$(_md)/bin/model_integrate
MODEL_INTEGRATE_BINARY?=$(_md)/bin.$(shell hostname)/model_integrate

# expected counts
$(EE_MAP_EXPECTED): $(EE_MAP_BINNED)
	$(_start)
	rm -rf $(EE_TMP_CONTIG_DIR) $(EE_LOG_CONTIG_DIR)
	$(_R) R/model_predict.r compute.expected.counts \
		   binary=$(MODEL_INTEGRATE_BINARY) \
		   tmp.dir=$(EE_TMP_CONTIG_DIR) \
		   log.dir=$(EE_LOG_CONTIG_DIR) \
		   wd=$(_md) \
	           model.prefix=$(EE_MAP_IN_MODEL_PREFIX) \
	           fends.ifn=$(EE_MAP_BINNED) \
	  	   model.ifn=$(EE_MAP_IN_MFN) \
		   scope=anchored model=std \
	           ofields.x=cluster_bin ofields.y=cluster_bin omit.x.zero=F omit.y.zero=F \
	  	   num.splits=$(EE_MAP_NUM_SPLIT) max.jobs.fn=$(MAX_JOBS_FN) req.mem=2000 dtype=$(DTYPE) \
		   Rcall="$(_Rcall)" \
	    	   ofn=$@
	$(_end)
ee_exp: $(EE_MAP_EXPECTED)

########################################################################################################################
# gene / contig / anchor map
########################################################################################################################

$(EE_MATRIX): $(EE_MAP_OBSERVED) $(EE_MAP_EXPECTED)
	$(_start)
	$(_md)/pl/ee_unite.pl \
		$(EE_MAP_OBSERVED) \
		$(EE_MAP_EXPECTED) \
		$(EE_MAP_BINNED).cluster \
		$@
	$(_end)
ee_map: $(EE_MATRIX)
