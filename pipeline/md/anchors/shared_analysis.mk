
#####################################################################################################
# generate tables with all genome tuples
#####################################################################################################

GS_REF_DONE?=$(GS_AAI_SHARED_DIR)/.done_ref
$(GS_REF_DONE):
	$(call _start,$(GS_AAI_SHARED_DIR))
	perl $(_md)/pl/aai_shared.pl \
		$(GS_REF_GENE_TABLE) \
		$(GS_REF_GENE_CLUSTER_TABLE) \
		$(GS_REF_GENE_FIELD) \
		$(GS_REF_GENOME_FIELD) \
		$(GS_REF_GENOME_TABLE) \
		$(GS_REF_AAI_TABLE) \
		$(GS_MIN_TUPLE_SIZE) \
		$(GS_MAX_TUPLE_SIZE) \
		$(GS_REF_PREFIX)
	$(_end_touch)
gs_ref: $(GS_REF_DONE)

GS_ANCHOR_DONE?=$(GS_AAI_SHARED_DIR)/.done_anchor
$(GS_ANCHOR_DONE):
	$(call _start,$(GS_AAI_SHARED_DIR))
	perl $(_md)/pl/aai_shared.pl \
		$(GS_ANCHOR_GENE_TABLE) \
		$(GS_ANCHOR_GENE_CLUSTER_TABLE) \
		$(GS_ANCHOR_GENE_FIELD) \
		$(GS_ANCHOR_GENOME_FIELD) \
		$(GS_ANCHOR_GENOME_TABLE) \
		$(GS_ANCHOR_AAI_TABLE) \
		$(GS_MIN_TUPLE_SIZE) \
		$(GS_MAX_TUPLE_SIZE) \
		$(GS_ANCHOR_PREFIX)
	$(_end_touch)
gs_anchor: $(GS_ANCHOR_DONE)

#####################################################################################################
# bin by AAI
#####################################################################################################

GS_REF_BIN_DONE?=$(GS_ANCHOR_BIN_DIR)/.done_ref_bin
$(GS_REF_BIN_DONE): $(GS_REF_DONE)
	$(call _start,$(GS_REF_BIN_DIR))
	$(_R) R/gs_analysis.r bin.gs.by.aai \
		ifn.prefix=$(GS_REF_PREFIX) \
		min.k=$(GS_MIN_TUPLE_SIZE) \
		max.k=$(GS_MAX_TUPLE_SIZE) \
		aai.breaks=$(GS_AAI_BREAKS) \
		odir=$(GS_REF_BIN_DIR)
	$(_end_touch)
gs_ref_bin: $(GS_REF_BIN_DONE)

GS_ANCHOR_BIN_DONE?=$(GS_ANCHOR_BIN_DIR)/.done_anchor_bin
$(GS_ANCHOR_BIN_DONE): $(GS_ANCHOR_DONE)
	$(call _start,$(GS_ANCHOR_BIN_DIR))
	$(_R) R/gs_analysis.r bin.gs.by.aai \
		ifn.prefix=$(GS_ANCHOR_PREFIX) \
		min.k=$(GS_MIN_TUPLE_SIZE) \
		max.k=$(GS_MAX_TUPLE_SIZE) \
		aai.breaks=$(GS_AAI_BREAKS) \
		odir=$(GS_ANCHOR_BIN_DIR)
	$(_end_touch)
gs_anchor_bin: $(GS_ANCHOR_BIN_DONE)

gs_all: $(GS_ANCHOR_BIN_DONE) $(GS_REF_BIN_DONE)

#####################################################################################################
# plots
#####################################################################################################

gs_plot:
	$(_start)
	$(_R) R/gs_analysis.r plot.gs \
		ref.dir=$(GS_REF_BIN_DIR) \
		anchor.dir=$(GS_ANCHOR_BIN_DIR) \
		min.k=$(GS_MIN_TUPLE_SIZE) \
		max.k=$(GS_MAX_TUPLE_SIZE) \
		aai.breaks=$(GS_AAI_BREAKS) \
		fdir=$(CA_MAP_FDIR)/shared_analysis/$(GR_SAMPLE_ID)/$(GS_AAI_BREAKS_ID)
	$(_end)
