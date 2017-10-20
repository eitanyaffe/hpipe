# fends
CA_MAP_FENDS?=$(CA_MAP_DIR)/fends
CA_MAP_BINNED?=$(CA_MAP_DIR)/fends.binned

CA_MAP_BINNED_FINAL?=$(CA_MAP_DIR)/fends.final

# observed
CA_MAP_OBSERVED?=$(CA_MAP_DIR)/count.o

# expected
CA_MAP_CONTIG_EXPECTED?=$(CA_MAP_DIR)/contig.e

# relative enr
CA_MAP_EN?=$(CA_MAP_DIR)/enrich


########################################################################################################################
# prepare fends table
########################################################################################################################

$(CA_MAP_FENDS):
	$(call _start,$(CA_MAP_DIR))
	$(_md)/pl/ca_fends.pl \
		$(CA_MAP_IN_FENDS) \
		$(CONTIG_TABLE) \
		$(ANCHOR_TABLE) \
		$@
	$(_end)
assign: $(CA_MAP_FENDS)

$(CA_MAP_BINNED): $(CA_MAP_FENDS)
	$(_start)
	$(_md)/pl/bin_fields.pl \
		$(CA_MAP_FENDS) \
		$(CA_MAP_BINNED) \
		contig anchor
	$(_end)
ca_binned: $(CA_MAP_BINNED)

$(CA_MAP_BINNED_FINAL): $(CA_MAP_BINNED)
	$(_start)
	$(_md)/pl/ca_fix_bins.pl $(CA_MAP_BINNED) $(CA_MAP_BINNED_FINAL)
	$(_end)
ca_fends: $(CA_MAP_BINNED_FINAL)

########################################################################################################################
# observed
########################################################################################################################

$(CA_MAP_OBSERVED): $(CA_MAP_BINNED)
	$(_start)
	$(_md)/pl/ca_observed.pl $(CA_MAP_BINNED) $(CA_MAP_IN_MAT) $(CA_CONTIG_NBINS) $@
	$(_end)
ca_obs: $(CA_MAP_OBSERVED)

########################################################################################################################
# expected contigs
########################################################################################################################

CA_QSUB_DIR?=$(QSUB_ANCHOR_DIR)/ca_map

TMP_CONTIG_DIR?=$(CA_QSUB_DIR)/contig_tmp
LOG_CONTIG_DIR?=$(CA_QSUB_DIR)/contig_log

MODEL_INTEGRATE_BINARY?=$(_md)/bin/model_integrate

# expected counts
$(CA_MAP_CONTIG_EXPECTED): $(CA_MAP_BINNED_FINAL)
	$(_start)
	rm -rf $(TMP_CONTIG_DIR) $(LOG_CONTIG_DIR)
	$(_R) R/model_predict.r compute.expected.counts \
		   binary=$(MODEL_INTEGRATE_BINARY) \
		   tmp.dir=$(TMP_CONTIG_DIR) \
		   log.dir=$(LOG_CONTIG_DIR) \
		   wd=$(_md) \
	           ifn.prefix=$(CA_MAP_IN_MODEL_PREFIX) \
		   cis.decay.table.dir=none \
	           fends.ifn=$(CA_MAP_BINNED_FINAL) \
	  	   model.ifn=$(CA_MAP_IN_MFN) \
		   filter=inter inter.field=contig_bin \
		   cis.threshold=0 compaction.fn=none cis.mode=none cis.binsize=0 \
	           ofields.x=contig_bin ofields.y=anchor_bin omit.x.zero=T omit.y.zero=T \
	  	   num.splits=$(CA_MAP_NUM_SPLIT) max.jobs.fn=$(MAX_JOBS_FN) req.mem=2000 dtype=$(DTYPE) \
		   Rcall="$(_Rcall)" \
	    	   ofn=$@
	$(_end)
ca_exp: $(CA_MAP_CONTIG_EXPECTED)

########################################################################################################################
# gene / contig / anchor map
########################################################################################################################

$(CA_MATRIX): $(CA_MAP_OBSERVED) $(CA_MAP_GENE_EXPECTED) $(CA_MAP_CONTIG_EXPECTED)
	$(_start)
	$(_md)/pl/ca_unite.pl \
		$(CA_MAP_OBSERVED) \
		$(CA_MAP_CONTIG_EXPECTED) \
		$(CA_MAP_BINNED).contig \
		$(CA_MAP_BINNED).anchor \
		$@
	$(_end)
ca_map: $(CA_MATRIX)

########################################################################################################################
# inter anchor matrix
########################################################################################################################

$(CA_INTER_ANCHOR_MATRIX): $(CA_MAP_OBSERVED) $(CA_MAP_GENE_EXPECTED) $(CA_MAP_CONTIG_EXPECTED)
	$(_start)
	$(_md)/pl/ca_inter_anchor_matrix.pl \
		$(CA_MAP_OBSERVED) \
		$(CA_MAP_CONTIG_EXPECTED) \
		$(CA_MAP_BINNED).contig \
		$(CA_MAP_BINNED).anchor \
		$(ANCHOR_TABLE) \
		$(CA_MIN_ANCHOR_CONTIGS) \
		$(CA_CONTIG_COVERAGE) \
		$@
	$(_end)
ca_inter_matrix: $(CA_INTER_ANCHOR_MATRIX)

########################################################################################################################
# genome cloud
########################################################################################################################

$(CA_ANCHOR_CONTIGS): $(CA_MATRIX)
	$(_start)
	$(_R) R/ca_analysis.r anchor.contigs \
		ifn=$(CA_MATRIX) \
		min.contacts=$(CA_MIN_CONTACTS) \
		min.enrichment=$(CA_MIN_ENRICHMENT) \
		min.anchor.contigs=$(CA_MIN_ANCHOR_CONTIGS) \
		min.contig.coverage=$(CA_CONTIG_COVERAGE) \
		fdr=$(CA_ASSIGN_FDR) \
		ofn=$@
	$(_end)

$(CA_ANCHOR_GENES): $(CA_ANCHOR_CONTIGS)
	$(_start)
	perl $(_md)/pl/anchor_genes.pl \
		$(CA_ANCHOR_CONTIGS) \
		$(GENE_TABLE) \
		$@
	$(_end)


make_ca: $(CA_ANCHOR_GENES) $(CA_INTER_ANCHOR_MATRIX)

########################################################################################################################
# plots
########################################################################################################################

plot_ca:
	$(_start)
	$(_R) R/ca_analysis.r plot.ca \
		ifn=$(CA_MATRIX) \
		fdir=$(CA_MAP_FDIR)/distribution \
		min.contacts=$(CA_MIN_CONTACTS) \
		min.enrichment=$(CA_MIN_ENRICHMENT) \
		min.anchor.contigs=$(CA_MIN_ANCHOR_CONTIGS) \
		min.contig.coverage=$(CA_CONTIG_COVERAGE) \
		fdr=$(CA_ASSIGN_FDR)
	$(_end)

# single plot with all anchor scatters
plot_ca_summary:
	$(_start)
	$(_R) R/ca_summary.r plot.ca.summary \
		ifn=$(CA_MATRIX) \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		min.contacts=$(CA_MIN_CONTACTS) \
		min.enrichment=$(CA_MIN_ENRICHMENT) \
		min.anchor.contigs=$(CA_MIN_ANCHOR_CONTIGS) \
		min.contig.coverage=$(CA_CONTIG_COVERAGE) \
		fdr=$(CA_ASSIGN_FDR) \
		fdir=$(CA_MAP_FDIR)/rotated_anchor_maps
	$(_end)

plot_inter_anchor:
	$(_R) R/ca_analysis.r plot.inter.anchor.matrix \
		ifn=$(CA_INTER_ANCHOR_MATRIX) \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		min.contacts=$(CA_MIN_CONTACTS) \
		min.enrichment=$(CA_MIN_ENRICHMENT) \
		fdir=$(ANCHOR_FIGURE_DIR)/inter_anchor

# plot anchor coverage
plot_cc_coverage: $(CA_ANCHOR_CONTIGS)
	$(_R) R/ca_analysis.r plot.cc.coverage \
		ifn=$(CA_ANCHOR_CONTIGS) \
		ifn.coverage=$(COVERAGE_TABLE) \
		fdir=$(CA_MAP_FDIR)/coverage

plot_abundance_distribution:
	$(_R) $(_md)/R/plot_abundance.r plot.abundance \
		ifn.coverage=$(COVERAGE_TABLE) \
		ifn.seeds=$(INITIAL_ANCHOR_TABLE) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.gc=$(CONTIG_GC_TABLE) \
		fdir=$(CA_MAP_FDIR)/assembly_abundance

plot_anchor_abundance:
	$(_R) $(_md)/R/plot_abundance.r plot.anchor.abundance \
		ifn=$(ANCHOR_TABLE) \
		ifn.coverage=$(COVERAGE_TABLE) \
		fdir=$(CA_MAP_FDIR)/anchor_abundance

plot_shared:
	$(_R) $(_md)/R/plot_shared.r plot.shared \
		ifn=$(CA_ANCHOR_CONTIGS) \
		ifn.contigs=$(CONTIG_TABLE) \
		fdir=$(CA_MAP_FDIR)/shared

ca_stats:
	$(_R) $(_md)/R/anchor_stats.r anchor.stats \
		ifn=$(CA_ANCHOR_CONTIGS) \
		ifn.contigs=$(CONTIG_TABLE) \
		odir=$(CA_MAP_FDIR)

make_ca_plots: ca_stats plot_ca plot_ca_summary plot_inter_anchor plot_cc_coverage plot_abundance_distribution plot_anchor_abundance plot_shared

########################################################################################################################
# [deprecated] element-anchor matrix, depends on repsonse clusters
########################################################################################################################

$(EA_MATRIX): $(CA_MAP_OBSERVED) $(CA_MAP_GENE_EXPECTED) $(CA_MAP_CONTIG_EXPECTED)
	$(_start)
	$(_md)/pl/ea_unite.pl \
		$(CA_MAP_OBSERVED) \
		$(CA_MAP_CONTIG_EXPECTED) \
		$(CA_MAP_BINNED).contig \
		$(CA_MAP_BINNED).anchor \
		$(RESPONSE_CONTIG_MEAN_CLUSTER) \
		$@
	$(_end)
ea_map: $(EA_MATRIX)

