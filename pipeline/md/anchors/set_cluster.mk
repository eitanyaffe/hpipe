SET_CLUSTER_MATRIX_DIR?=$(SET_CLUSTER_DIR)/map.dir
SET_CLUSTER_PROJECT_DIR?=$(SET_CLUSTER_DIR)/project.dir
SET_CLUSTER_MATRIX_TABLE?=$(SET_CLUSTER_DIR)/summary.table

# skip if set is equal to this value
SET_CLUSTER_NULL_VALUE?=0

####################################################################################
# one file per each set-set pair
####################################################################################

SET_CLUSTER_MATRIX_FILES_DONE?=$(SET_CLUSTER_DIR)/.done_split
$(SET_CLUSTER_MATRIX_FILES_DONE):
	$(call _start,$(SET_CLUSTER_MATRIX_DIR))
	$(_md)/pl/split_gene_map.pl \
		$(SET_CLUSTER_INPUT_MATRIX) \
		$(SET_CLUSTER_INPUT_FIELD) \
		$(SET_CLUSTER_NULL_VALUE) \
		$(SET_CLUSTER_INPUT_MATRIX) \
		$(SET_CLUSTER_INPUT_FIELD) \
		$(SET_CLUSTER_NULL_VALUE) \
		$(SET_CLUSTER_INPUT_MAP) \
		$(SET_CLUSTER_MATRIX_DIR)
	$(_end_touch)

SET_CLUSTER_MATRIX_PROJECT_DONE?=$(SET_CLUSTER_DIR)/.done_project
PROJECT_SCRIPT?=$(_md)/pl/map_to_functions.pl
$(SET_CLUSTER_MATRIX_PROJECT_DONE): $(SET_CLUSTER_MATRIX_FILES_DONE)
	$(call _start,$(SET_CLUSTER_PROJECT_DIR))
	$(_R) R/set_matrix.r distrib.top.hit \
		script=$(PROJECT_SCRIPT) \
		idir=$(SET_CLUSTER_MATRIX_DIR) \
		odir=$(SET_CLUSTER_PROJECT_DIR) \
		qsub.dir=$(SET_CLUSTER_QSUB_DIR) \
		batch.max.jobs=$(MAX_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		self=$(SET_SELF) \
		jobname=map2single
	$(_end_touch)

####################################################################################
# compute matrix and cluster it
####################################################################################

SET_CLUSTER_MATRIX_TABLE_DONE?=$(SET_CLUSTER_DIR)/.done_table
$(SET_CLUSTER_MATRIX_TABLE_DONE): $(SET_CLUSTER_MATRIX_PROJECT_DONE)
	$(_start)
	$(_R) R/set_matrix.r compute.matrix.single \
		mat.dir=$(SET_CLUSTER_PROJECT_DIR) \
		ifn=$(SET_CLUSTER_INPUT_MATRIX) \
		field=$(SET_CLUSTER_INPUT_FIELD) \
		ifn.genes=$(SET_GENES_TABLE) \
		identity.threshold=$(MEAN_IDENTITY_THRESHOLD) \
		coverage.threshold=$(COVERAGE_THRESHOLD) \
		min.gene.count=$(MIN_GENE_COUNT) \
		include.shared=$(INCLUDE_SHARED) \
		ofn=$(SET_CLUSTER_MATRIX_TABLE)
	$(_end_touch)

SET_CLUSTER_TABLE_DONE?=$(SET_CLUSTER_DIR)/.done_cluster
$(SET_CLUSTER_TABLE_DONE): $(SET_CLUSTER_MATRIX_TABLE_DONE)
	$(_start)
	$(_R) R/set_matrix.r cluster.sets \
		mat.ifn=$(SET_CLUSTER_MATRIX_TABLE) \
		prefix=$(SET_CLUSTER_PREFIX) \
		min.identity=$(MEAN_IDENTITY_THRESHOLD) \
		ofn=$(SET_CLUSTER_TABLE)
	$(_end_touch)

####################################################################################
# make
####################################################################################

make_set_clusters: $(SET_CLUSTER_TABLE_DONE)

####################################################################################
# plot
####################################################################################

plot_set_clusters:
	$(_start)
	$(_R) R/set_matrix.r plot.matrix.single \
		mat.ifn=$(SET_CLUSTER_MATRIX_TABLE) \
		min.identity=$(MEAN_IDENTITY_THRESHOLD) \
		title=$(SET_TITLE) \
		order.ifn=$(SET_CLUSTER_TABLE) \
		fdir=$(SET_CLUSTER_FDIR)
	$(_end)

plot_set_cluster_full:
	$(_start)
	$(_R) R/set_matrix.r plot.matrix.full.single \
		mat.ifn=$(SET_CLUSTER_MATRIX_TABLE) \
		title=$(SET_TITLE) \
		order.ifn=$(SET_CLUSTER_TABLE) \
		fdir=$(SET_CLUSTER_FDIR)
	$(_end)

####################################################################################
# main
####################################################################################

.PHONY: plot_set_clusters plot_set_cluster_full
