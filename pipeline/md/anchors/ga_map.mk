# fends
GA_MAP_FENDS?=$(GA_MAP_DIR)/fends
GA_MAP_BINNED?=$(GA_MAP_DIR)/fends.binned

GA_MAP_BINNED_FINAL?=$(GA_MAP_DIR)/fends.final

# observed
GA_MAP_OBSERVED?=$(GA_MAP_DIR)/count.o

# expected
GA_MAP_GENE_EXPECTED?=$(GA_MAP_DIR)/gene.e
GA_MAP_CONTIG_EXPECTED?=$(GA_MAP_DIR)/contig.e

########################################################################################################################
# prepare fends table
########################################################################################################################

$(GA_MAP_FENDS):
	$(call _start,$(GA_MAP_DIR))
	$(_md)/pl/assign_anchor_to_fends.pl \
		$(GA_MAP_IN_FENDS) \
		$(CONTIG_TABLE) \
		$(ANCHOR_TABLE) \
		$(GENE_TABLE) \
		$@
	$(_end)
assign: $(GA_MAP_FENDS)

$(GA_MAP_BINNED): $(GA_MAP_FENDS)
	$(_start)
	$(_md)/pl/bin_fields.pl \
		$(GA_MAP_FENDS) \
		$(GA_MAP_BINNED) \
		contig gene anchor
	$(_end)
ga_binned: $(GA_MAP_BINNED)

$(GA_MAP_BINNED_FINAL): $(GA_MAP_BINNED)
	$(_start)
	$(_md)/pl/ga_fix_bins.pl $(GA_MAP_BINNED) $(GA_MAP_BINNED_FINAL)
	$(_end)
ga_fends: $(GA_MAP_BINNED_FINAL)

########################################################################################################################
# observed
########################################################################################################################

$(GA_MAP_OBSERVED): $(GA_MAP_BINNED)
	$(_start)
	$(_md)/pl/ga_observed.pl $(GA_MAP_BINNED) $(GA_MAP_IN_MAT) $(GA_CONTIG_NBINS) $@
	$(_end)
obs: $(GA_MAP_OBSERVED)

########################################################################################################################
# expected contigs
########################################################################################################################

GA_QSUB_DIR?=$(QSUB_DIR)/datasets/$(DATASET)/ga_map

TMP_CONTIG_DIR?=$(GA_QSUB_DIR)/contig_tmp
LOG_CONTIG_DIR?=$(GA_QSUB_DIR)/contig_log

MODEL_INTEGRATE_BINARY?=$(_md)/bin/model_integrate

# expected counts
$(GA_MAP_CONTIG_EXPECTED): $(GA_MAP_BINNED_FINAL)
	$(_start)
	rm -rf $(TMP_CONTIG_DIR) $(LOG_CONTIG_DIR)
	$(_R) R/model_predict.r compute.expected.counts \
		   binary=$(MODEL_INTEGRATE_BINARY) \
		   tmp.dir=$(TMP_CONTIG_DIR) \
		   log.dir=$(LOG_CONTIG_DIR) \
		   wd=$(_md) \
	           ifn.prefix=$(GA_MAP_IN_MODEL_PREFIX) \
		   cis.decay.table.dir=none \
	           fends.ifn=$(GA_MAP_BINNED_FINAL) \
	  	   model.ifn=$(GA_MAP_IN_MFN) \
		   filter=inter inter.field=contig_bin \
		   cis.threshold=0 compaction.fn=none cis.mode=none cis.binsize=0 \
	           ofields.x=contig_bin ofields.y=anchor_bin omit.x.zero=T omit.y.zero=T \
	  	   num.splits=$(GA_MAP_NUM_SPLIT) max.jobs.fn=$(MAX_JOBS_FN) req.mem=2000 dtype=$(DTYPE) \
		   Rcall="$(_Rcall)" \
	    	   ofn=$@
	$(_end)
exp_contigs: $(GA_MAP_CONTIG_EXPECTED)

########################################################################################################################
# expected genes
########################################################################################################################

TMP_GENE_DIR?=$(GA_QSUB_DIR)/gene_tmp
LOG_GENE_DIR?=$(GA_QSUB_DIR)/gene_log

# expected counts
$(GA_MAP_GENE_EXPECTED): $(GA_MAP_BINNED_FINAL)
	$(_start)
	rm -rf $(TMP_GENE_DIR) $(LOG_GENE_DIR)
	$(_R) R/model_predict.r compute.expected.counts \
		   binary=$(MODEL_INTEGRATE_BINARY) \
		   tmp.dir=$(TMP_GENE_DIR) \
		   log.dir=$(LOG_GENE_DIR) \
		   wd=$(_md) \
	           ifn.prefix=$(GA_MAP_IN_MODEL_PREFIX) \
		   cis.decay.table.dir=none \
	           fends.ifn=$(GA_MAP_BINNED_FINAL) \
	  	   model.ifn=$(GA_MAP_IN_MFN) \
		   filter=inter inter.field=contig_bin \
		   cis.threshold=0 compaction.fn=none cis.mode=none cis.binsize=0 \
	           ofields.x=gene_bin ofields.y=anchor_bin omit.x.zero=T omit.y.zero=T \
	  	   num.splits=$(GA_MAP_NUM_SPLIT) max.jobs.fn=$(MAX_JOBS_FN) req.mem=2000 dtype=$(DTYPE) \
		   Rcall="$(_Rcall)" \
	    	   ofn=$@
	$(_end)
exp_genes: $(GA_MAP_GENE_EXPECTED)

########################################################################################################################
# gene / contig / anchor map
########################################################################################################################

$(GA_MATRIX): $(GA_MAP_OBSERVED) $(GA_MAP_GENE_EXPECTED) $(GA_MAP_CONTIG_EXPECTED)
	$(_start)
	$(_md)/pl/ga_unite.pl \
		$(GA_MAP_OBSERVED) \
		$(GA_MAP_GENE_EXPECTED) \
		$(GA_MAP_CONTIG_EXPECTED) \
		$(GA_MAP_BINNED).gene \
		$(GA_MAP_BINNED).contig \
		$(GA_MAP_BINNED).anchor \
		$@
	$(_end)
ga_map: $(GA_MATRIX)

########################################################################################################################
# extending anchors
########################################################################################################################

GA_MAP_FDIR?=$(ANCHOR_FIGURE_DIR)/ga_map

# plot to display threshold and complete scatters
plot_ga:
	$(_start)
	$(_R) R/ga_analysis.r plot.ga.distrib \
		ifn=$(GA_MATRIX) \
		fdir=$(GA_MAP_FDIR)/distribution \
		min.contacts=$(GA_MIN_CONTACTS) \
		min.anchor.contigs=$(GA_MIN_ANCHOR_CONTIGS) \
		min.contig.coverage=$(GA_CONTIG_COVERAGE) \
		q.threshold=$(GA_Q_THRESHOLD)
	$(_end)

#		marked.contigs.ifn=$(MARKED_CONTIGS) \

plot_ga_matrix:
	$(_start)
	$(_R) R/gene_graph.r plot.ga.matrix \
		ifn.ga=$(GA_ANCHOR_GENES) \
		fdir=$(GA_MAP_FDIR)/matrix
	$(_end)

$(GA_ANCHOR_GENES): $(GA_MATRIX)
	$(_start)
	$(_R) R/ga_analysis.r anchor.genes \
		ifn=$(GA_MATRIX) \
		ifn.cgenes=$(GENE_CLUSTER_TABLE) \
		min.contacts=$(GA_MIN_CONTACTS) \
		min.anchor.contigs=$(GA_MIN_ANCHOR_CONTIGS) \
		min.contig.coverage=$(GA_CONTIG_COVERAGE) \
		q.threshold=$(GA_Q_THRESHOLD) \
		ofn=$@
	$(_end)
make_ga: $(GA_ANCHOR_GENES)

########################################################################################################################
# compute various genomic stats
########################################################################################################################

# info:
# - contig median size
# - average gc content
# - number of shared genes
# - relative abundance
# - contact enrichment
# - number of genes
# - number of contigs
ANCHOR_INFO_DONE?=$(GA_MAP_DIR)/.done_info
$(ANCHOR_INFO_DONE):
	$(_start)
	$(_R) R/ga_analysis.r anchor.info \
		ga.ifn=$(GA_ANCHOR_GENES) \
		gene.ifn=$(GENE_TABLE) \
		contig.ifn=$(CONTIG_TABLE) \
		gc.ifn=$(CONTIG_GC_TABLE) \
		coverage.ifn=$(COVERAGE_TABLE) \
		ofn.genes=$(ANCHOR_INFO_TABLE_GENES) \
		ofn.contigs=$(ANCHOR_INFO_TABLE_CONTIGS) \
		ofn.summary=$(ANCHOR_INFO_TABLE)
	$(_end_touch)
ainfo: $(ANCHOR_INFO_DONE)

########################################################################################################################
# plots
########################################################################################################################

plot_anchor_info: $(ANCHOR_INFO_DONE)
	$(_start)
	$(_R) R/ga_analysis.r plot.anchor.info \
		ifn=$(ANCHOR_INFO_TABLE) \
		ifn.genes=$(ANCHOR_INFO_TABLE_GENES) \
		ifn.contigs=$(ANCHOR_INFO_TABLE_CONTIGS) \
		fdir=$(GA_MAP_FDIR)/info
	$(_R) R/ga_analysis.r plot.anchor.info.summary \
		ifn=$(ANCHOR_INFO_TABLE) \
		fdir=$(GA_MAP_FDIR)/info
	$(_end)

plot_anchor_coverage:
	$(_start)
	$(_R) R/ga_analysis.r anchor.coverage \
		ifn=$(GA_ANCHOR_GENES) \
		ifn.cov=$(COVERAGE_TABLE) \
		fdir=$(GA_MAP_FDIR)/coverage
	$(_end)

########################################################################################################################
# gene graph
########################################################################################################################

plot_gene_graph:
	$(_start)
	$(_R) R/gene_graph.r plot.gene.graph \
		ifn.genes=$(GENE_TABLE) \
		ifn.ga=$(GA_ANCHOR_GENES) \
		fdir=$(GA_MAP_FDIR)
	$(_end)

make_ga_plots: plot_ga plot_anchor_coverage plot_anchor_info
