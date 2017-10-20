########################################################################################################################
# taxa analysis
########################################################################################################################

# clustering
TAXA_MATRIX?=$(TAXA_DIR)/matrix
TAXA_SET_ORDER?=$(TAXA_DIR)/anchor_order

TAXA_TREES_DONE?=$(TAXA_DIR)/.done_trees_taxa
$(TAXA_TREES_DONE):
	$(call _start,$(TAXA_DIR))
	$(call _time,$(TAXA_DIR),trees) \
		$(_md)/pl/compute_tax_trees.pl \
			$(TAXA_SET_GENES) \
			$(TAXA_SET_FIELD) \
			$(UNIREF_TOP) \
			$(UNIREF_TAX_LOOKUP) \
			$(NCBI_TAX_NODES) \
			$(NCBI_TAX_NAMES) \
			$(NCBI_TAX_MERGED) \
			$(NCBI_TAX_DELETED) \
			$(TAX_MIN_IDENTITY) \
			$(TAX_MIN_COVER) \
			$(CAG_MASK) \
			$(TAX_MASK) \
			$(SET_TAXA_SUMMARY) \
			$(SET_TAXA_TABLE) \
			$(SET_TAXA_TREES)
	$(_end_touch)
taxa_trees: $(TAXA_TREES_DONE)

# identify path leading to leaf with maximal coverage in terms of genes
TAXA_REP_PATH_DONE?=$(TAXA_DIR)/.done_rep_path
$(TAXA_REP_PATH_DONE): $(TAXA_TREES_DONE)
	$(_start)
	$(_R) R/taxa_path.r taxa.path \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.tree=$(SET_TAXA_TREES) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.gene.table=$(GENEBANK_TABLE) \
		ifn.summary=$(SET_TAXA_SUMMARY) \
		ofn.path=$(SET_TAXA_PATH) \
		ofn.species=$(SET_TAXA_REPS)
	$(_end_touch)
taxa_path: $(TAXA_REP_PATH_DONE)

# determine resolution
TAXA_RESOLVE_DONE?=$(TAXA_DIR)/.done_resolve
$(TAXA_RESOLVE_DONE): $(TAXA_REP_PATH_DONE)
	$(_start)
	$(_R) R/taxa_resolve.r taxa.resolve \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.path=$(SET_TAXA_PATH) \
		ifn.tree=$(SET_TAXA_TREES) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.summary=$(SET_TAXA_SUMMARY) \
		min.coverage=$(SET_TAXA_RESOLVE_MIN_COVERAGE) \
		min.ratio=$(SET_TAXA_RESOLVE_MIN_RATIO) \
		ofn=$(SET_TAXA_RESOLVE)
	$(_end_touch)
taxa_resolve: $(TAXA_RESOLVE_DONE)

TAXA_REP_LEVELS_DONE?=$(TAXA_DIR)/.done_levels
$(TAXA_REP_LEVELS_DONE): $(TAXA_REP_PATH_DONE)
	$(_start)
	$(_R) R/taxa_legend.r taxa.levels \
		ifn.species=$(SET_TAXA_REPS) \
		ifn.path=$(SET_TAXA_PATH) \
		levels=$(SET_TAX_LEVELS) \
		ofn.level.prefix=$(SET_TAX_LEVEL_PREFIX) \
		ofn.level.summary.prefix=$(SET_TAX_LEVEL_SUMMARY_PREFIX)
	$(_end_touch)
taxa_levels: $(TAXA_REP_LEVELS_DONE)

TAXA_REP_LEGEND_DONE?=$(TAXA_DIR)/.done_legend
$(TAXA_REP_LEGEND_DONE): $(TAXA_REP_LEVELS_DONE)
	$(_start)
	$(_R) R/taxa_legend.r taxa.legend \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.species=$(SET_TAXA_REPS) \
		ifn.level.prefix=$(SET_TAX_LEVEL_PREFIX) \
		ifn.level.summary.prefix=$(SET_TAX_LEVEL_SUMMARY_PREFIX) \
		letter.level=$(SET_TAX_LETTER_LEVEL) \
		letter.order=$(SET_TAX_LETTER_ORDER) \
		force.ids=$(SET_TAX_LEGEND_IDS) \
		force.colors=$(SET_TAX_LEGEND_COLORS) \
		color.level=$(SET_TAX_COLOR_LEVEL) \
		ofn.legend.letter=$(SET_TAX_LEGEND_LETTER) \
		ofn.legend.color=$(SET_TAX_LEGEND_COLOR) \
		ofn=$(SET_TAX_LEGEND)
	$(_end_touch)
taxa_legend: $(TAXA_REP_LEGEND_DONE)

TAXA_GENOMES_DONE?=$(TAXA_DIR)/.done_genomes
$(TAXA_GENOMES_DONE): $(TAXA_TREES_DONE)
	$(_start)
	$(_R) R/taxa_reps.r select.species \
		ifn=$(SET_TAXA_TREES) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.gene.table=$(GENEBANK_TABLE) \
		level=$(TAXA_REF_LEVEL) \
		ofn=$(TAXA_REF_GENOMES)
	$(_end_touch)
taxa_genomes: $(TAXA_GENOMES_DONE)

.PHONY: taxa_trees taxa_reps GO_append make_taxa

###################################################################################################
# subtree
###################################################################################################

# generate subtrees for a list of ids
TAXA_SUBTREES_DONE?=$(TAXA_SUBTREE_DIR)/.done_sub_trees
$(TAXA_SUBTREES_DONE):
	$(call _start,$(TAXA_SUBTREE_DIR))
	$(_md)/pl/subtree.pl \
		$(TAXA_SUBTREES_INPUT) \
		$(NCBI_TAX_NODES) \
		$(NCBI_TAX_NAMES) \
		$(NCBI_TAX_MERGED) \
		$(NCBI_TAX_DELETED) \
		$(TAXA_SUBTREES_TABLE)
	$(_end_touch)
taxa_subtrees: $(TAXA_SUBTREES_DONE)

###################################################################################################
# old code for ref representatives
###################################################################################################

# select rep taxa node by going down tree until ambigous split
TAXA_REP_DONE?=$(TAXA_DIR)/.done_rep
$(TAXA_REP_DONE): $(TAXA_TREES_DONE)
	$(_start)
	$(_R) R/taxa_reps.r select.reps \
		ifn=$(SET_TAXA_TREES) \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.summary=$(SET_TAXA_SUMMARY) \
		min.ratio.1=$(TAXA_SELECT_RATIO_1) \
		min.ratio.12=$(TAXA_SELECT_RATIO_12) \
		skip.ids=$(TAX_ID_SKIP) \
		skip.names=$(TAX_NAME_SKIP) \
		ofn=$(SET_TAXA_REPS)
	$(_end_touch)
taxa_rep: $(TAXA_REP_DONE)
