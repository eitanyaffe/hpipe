
# skip if set is equal to this value
SET_COMPARE_NULL_VALUE1?=0
SET_COMPARE_NULL_VALUE2?=0

####################################################################################
# one file per each ref-anchor pair
####################################################################################

SET_COMPARE_GENE_MATRIX_DONE?=$(SET_COMPARE_DIR)/.done_split
$(SET_COMPARE_GENE_MATRIX_DONE): $(REF_CLUSTER_TABLE) $(ANCHOR_CLUSTER_TABLE)
	$(call _start,$(SET_COMPARE_MATRIX_DIR))
	$(_md)/pl/split_gene_map.pl \
		$(REF_SET_INPUT1) \
		$(REF_FIELD_INPUT1) \
		$(SET_COMPARE_NULL_VALUE1) \
		$(REF_SET_INPUT2) \
		$(REF_FIELD_INPUT2) \
		$(SET_COMPARE_NULL_VALUE2) \
		$(REF_MAP_INPUT) \
		$(SET_COMPARE_MATRIX_DIR)
	$(_end_touch)

SET_COMPARE_GENE_PROJECT_DONE?=$(SET_COMPARE_DIR)/.done_project
PROJECT_SCRIPT?=$(_md)/pl/map_to_functions.pl
$(SET_COMPARE_GENE_PROJECT_DONE): $(SET_COMPARE_GENE_MATRIX_DONE)
	$(call _start,$(SET_COMPARE_PROJECT_DIR))
	$(_R) R/set_matrix.r distrib.top.hit \
		script=$(PROJECT_SCRIPT) \
		idir=$(SET_COMPARE_MATRIX_DIR) \
		odir=$(SET_COMPARE_PROJECT_DIR) \
		qsub.dir=$(SET_COMPARE_QSUB) \
		batch.max.jobs=$(MAX_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		self=F \
		jobname=map2single
	$(_end_touch)

####################################################################################
# compute matrix
####################################################################################

$(SET_COMPARE_MATRIX_TABLE): $(SET_COMPARE_GENE_PROJECT_DONE)
	$(_start)
	$(_R) R/set_matrix.r compute.matrix \
		mat.dir=$(SET_COMPARE_PROJECT_DIR) \
		ifn1=$(REF_SET_INPUT1) \
		ifn.genes1=$(REF_SET_INPUT_GENES1) \
		field1=$(REF_FIELD_INPUT1) \
		ifn2=$(REF_SET_INPUT2) \
		ifn.genes2=$(REF_SET_INPUT_GENES2) \
		field2=$(REF_FIELD_INPUT2) \
		identity.threshold=$(MEAN_IDENTITY_THRESHOLD) \
		coverage.threshold=$(COVERAGE_THRESHOLD) \
		min.gene.count=$(MIN_GENE_COUNT) \
		include.shared=$(INCLUDE_SHARED) \
		ofn=$@
	$(_end)
sc_matrix: $(SET_COMPARE_MATRIX_TABLE)

####################################################################################
# project to all
####################################################################################

GENE_PROJECT_FWD_DONE?=$(SET_COMPARE_DIR)/.done_fwd
$(GENE_PROJECT_FWD_DONE):
	$(_start)
	perl $(_md)/pl/project_gene.pl \
		$(REF_SET_INPUT1) \
		$(REF_FIELD_INPUT1) \
		$(REF_MAP_INPUT) \
		gene1 gene2 \
		$(GENE_PROJECT_FWD)
	$(_end_touch)

GENE_PROJECT_BCK_DONE?=$(SET_COMPARE_DIR)/.done_bck
$(GENE_PROJECT_BCK_DONE): $(GENE_PROJECT_MAP)
	$(_start)
	perl $(_md)/pl/project_gene.pl \
		$(REF_SET_INPUT2) \
		$(REF_FIELD_INPUT2) \
		$(REF_MAP_INPUT) \
		gene2 gene1 \
		$(GENE_PROJECT_BCK)
	$(_end_touch)

project_all: $(GENE_PROJECT_FWD_DONE) $(GENE_PROJECT_BCK_DONE)

plot_basic_project: $(GENE_PROJECT_FWD_DONE) $(GENE_PROJECT_BCK_DONE)
	$(_start)
	$(_R) R/set_matrix.r plot.gene.project \
		order.ifn1=$(REF_CLUSTER_TABLE1) \
		order.ifn2=$(REF_CLUSTER_TABLE2) \
		fwd.source.ifn=$(GENE_PROJECT_FWD_SOURCE) \
		bck.source.ifn=$(GENE_PROJECT_BCK_SOURCE) \
		fdir=$(ANCHOR_FIGURE_DIR)/basic_gene_projection
	$(_end)

####################################################################################
# set best match
####################################################################################

SET_MATCH_DONE?=$(SET_COMPARE_DIR)/.done_match
$(SET_MATCH_DONE): $(SET_COMPARE_MATRIX_TABLE)
	$(_start)
	$(_R) R/set_matrix.r compute.best.match \
		ifn=$(SET_COMPARE_MATRIX_TABLE) \
		ofn1=$(SET_MATCH_TABLE1) \
		ofn2=$(SET_MATCH_TABLE2)
	$(_end_touch)

ref_matrix: $(SET_MATCH_DONE)

plot_match: $(SET_MATCH_DONE)
	$(_start)
	$(_R) R/set_matrix.r plot.match \
		ifn1=$(SET_MATCH_TABLE1) \
		ifn2=$(SET_MATCH_TABLE2) \
		fdir=$(SET_COMPARE_FDIR)
	$(_end)

####################################################################################
# project to best match.
# when projecting P back we include all R_i such that P = best_match(R_i)
####################################################################################

SET_PROJECT_DONE?=$(SET_COMPARE_DIR)/.done_gene_project
$(SET_PROJECT_DONE): $(SET_MATCH_DONE) $(GENE_PROJECT_FWD_DONE) $(GENE_PROJECT_BCK_DONE) $(SET_COMPARE_GENE_PROJECT_DONE)
	$(call _start,$(SET_TMP_DIR))
	$(_R) R/set_matrix.r project.genes \
		match.ifn1=$(SET_MATCH_TABLE1) \
		match.ifn2=$(SET_MATCH_TABLE2) \
		all.ifn1=$(GENE_PROJECT_FWD) \
		all.ifn2=$(GENE_PROJECT_BCK) \
		tmp.dir=$(SET_TMP_DIR) \
		collapse.script=$(_md)/pl/merge_collapse_maps.pl \
		mat.dir=$(SET_COMPARE_PROJECT_DIR) \
		ifn1=$(REF_SET_INPUT1) \
		ifn.genes1=$(REF_SET_INPUT_GENES1) \
		field1=$(REF_FIELD_INPUT1) \
		ifn2=$(REF_SET_INPUT2) \
		ifn.genes2=$(REF_SET_INPUT_GENES2) \
		field2=$(REF_FIELD_INPUT2) \
		ofn1=$(SET_FWD_TABLE) \
		ofn2=$(SET_BCK_TABLE)
	$(_end_touch)
set_project: $(SET_PROJECT_DONE)

####################################################################################
# plots
####################################################################################

plot_ref_matrix: $(SET_COMPARE_MATRIX_TABLE)
	$(_start)
	$(_R) R/set_matrix.r plot.matrix \
		mat.ifn=$(SET_COMPARE_MATRIX_TABLE) \
		min.identity=$(MEAN_IDENTITY_THRESHOLD) \
		order.ifn1=$(REF_CLUSTER_TABLE1) \
		order.ifn2=$(REF_CLUSTER_TABLE2) \
		title1=$(REF_SET_TITLE1) \
		title2=$(REF_SET_TITLE2) \
		fdir=$(SET_COMPARE_FDIR)
	$(_R) R/set_matrix.r plot.matrix.full \
		mat.ifn=$(SET_COMPARE_MATRIX_TABLE) \
		order.ifn1=$(REF_CLUSTER_TABLE1) \
		order.ifn2=$(REF_CLUSTER_TABLE2) \
		title1=$(REF_SET_TITLE1) \
		title2=$(REF_SET_TITLE2) \
		fdir=$(SET_COMPARE_FDIR)
	$(_end)

plot_project: $(SET_PROJECT_DONE)
	$(_start)
	$(_R) R/set_matrix.r plot.project \
		ifn.classify.contigs=$(CLASSIFY_CONTIG_TABLE) \
		ifn.classify.genes=$(CLASSIFY_GENE_TABLE) \
		ifn.classify.ref.genes=$(CLASSIFY_REF_GENE_TABLE) \
		ifn.fwd=$(SET_FWD_TABLE) \
		ifn.bck=$(SET_BCK_TABLE) \
		order.ifn.fwd=$(REF_CLUSTER_TABLE1) \
		order.ifn.bck=$(REF_CLUSTER_TABLE2) \
		fdir=$(ANCHOR_FIGURE_DIR)/ref_verify

# !!! plot_match requires igraph, so left out  
make_ref_plots: plot_ref_matrix plot_project

