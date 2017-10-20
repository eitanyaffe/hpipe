
copy_taxa_table:
	mkdir -p $(TAXA_FDIR)
	cp $(SC_SUMMARY_UNIQUE) $(SET_TAXA_REPS) $(TAXA_FDIR)

plot_taxa_uniref:
	$(_R) R/plot_taxa.r plot.uniref.summary \
		ifn.summary=$(SET_TAXA_SUMMARY) \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		fdir=$(TAXA_FDIR)/1_uniref

plot_taxa_path:
	$(_R) R/taxa_path_plot.r plot.path \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.path=$(SET_TAXA_PATH) \
		ifn.species=$(SET_TAXA_REPS) \
		fdir=$(TAXA_FDIR)/2.1_tree_path

# anchor trees
plot_taxa_trees:
	$(_R) R/plot_taxa_trees.r plot.anchor.trees \
		ifn.tree=$(SET_TAXA_TREES) \
		ifn.rep=$(SET_TAXA_REPS) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.res=$(SET_TAXA_RESOLVE) \
		ifn.summary=$(SET_TAXA_SUMMARY) \
		min.gene.count=$(SET_TAXA_TREE_PLOT_MIN_GENES) \
		fdir=$(TAXA_FDIR)/2.2_tree_whole

# anchor trees
plot_taxa_limited_trees:
	$(_R) R/plot_taxa_trees.r plot.anchor.trees.limited \
		ifn.tree=$(SET_TAXA_TREES) \
		ifn.rep=$(SET_TAXA_REPS) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.res=$(SET_TAXA_RESOLVE) \
		ifn.summary=$(SET_TAXA_SUMMARY) \
		min.gene.count=$(SET_TAXA_TREE_PLOT_MIN_GENES) \
		fdir=$(TAXA_FDIR)/2.3_tree_limited

plot_taxa_legend:
	$(_R) R/taxa_legend_plot.r plot.taxa.legend.color \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.legend.color=$(SET_TAX_LEGEND_COLOR) \
		fdir=$(TAXA_FDIR)/3_taxa_legend
	$(_R) R/taxa_legend_plot.r plot.taxa.legend.letter \
		ifn.taxa=$(SET_TAXA_TABLE) \
		ifn.legend.letter=$(SET_TAX_LEGEND_LETTER) \
		ifn.legend.color=$(SET_TAX_LEGEND_COLOR) \
		fdir=$(TAXA_FDIR)/3_taxa_legend

plot_taxa_rep_detail:
	$(_R) R/sc_map_plot.r plot.details \
		ifn.pairs=$(SC_SUMMARY_UNIQUE) \
		fdir=$(TAXA_FDIR)/4_rep_detail

# plot anchor/ref mapping distrib
plot_taxa_rep:
	$(_R) R/sc_map_plot.r plot.mapping \
		ifn.pairs=$(SC_SUMMARY_UNIQUE) \
		length=$(SC_FRAGMENT_LENGTH) \
		anchor2ref.dir=$(SC_ANCHOR2REF_DIR) \
		ref2anchor.dir=$(SC_REF2ANCHOR_DIR) \
		fdir=$(TAXA_FDIR)/4.1_rep

plot_taxa_rep_analysis:
	$(_R) R/sc_map_plot.r plot.rep.analysis \
		ifn.rep=$(SC_SUMMARY_UNIQUE) \
		ifn.resolve=$(SET_TAXA_RESOLVE) \
		checkm.ifn=$(CHECKM_QA) \
		min.completeness=$(SET_TAXA_MIN_COMPLETENESS) \
		fdir=$(TAXA_FDIR)/4.2_rep_analysis

# rep trees
plot_taxa_single_tree:
	$(_R) R/plot_taxa_reps.r plot.rep.tree \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.rep=$(SET_TAXA_REPS) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		trim.levels=$(TRIM_TREE_LEVELS) \
		trim.names=$(TRIM_TREE_NAMES) \
		trim.grep.names=$(TRIM_TREE_GREP_NAMES) \
		collapse.depth.ids=$(COLLAPSE_DEPTH_IDS) \
		level.names=$(LEVEL_NAMES) \
		level.cols=$(LEVEL_COLORS) \
		branch.ids=$(BRANCH_IDS) \
		branch.color.indices=$(BRANCH_COLOR_INDICES) \
		no.branch.color=$(NO_BRANCH_COLOR) \
		no.level.color=$(NO_LEVEL_COLOR) \
		ofn.legend=$(SET_TAXA_REP_LEGEND) \
		fdir=$(TAXA_FDIR)/5_rep_summary_tree

plot_taxa_summary_basic:
	$(_R) R/anchor_plot.r plot.complete.basic \
		contig.ifn=$(CONTIG_TABLE) \
		ca.ifn=$(CA_ANCHOR_CONTIGS) \
		coverage.ifn=$(ANCHOR_COVERAGE_TABLE) \
		gc.ifn=$(ANCHOR_GC_TABLE) \
		mat.ifn=$(ANCHOR_MATRIX_TABLE) \
		order.ifn=$(ANCHOR_CLUSTER_TABLE) \
		info.ifn=$(ANCHOR_SIZE_STATS) \
		group.ifn=$(GROUP_ANCHOR_TABLE) \
		taxa.ifn=$(SET_TAX_LEGEND) \
		checkm.ifn=$(CHECKM_QA) \
		min.identity=$(MEAN_IDENTITY_THRESHOLD) \
		fdir=$(TAXA_FDIR)/6_combined

plot_taxa_summary_taxa:
	$(_R) R/anchor_plot.r plot.complete.taxa \
		order.ifn=$(ANCHOR_CLUSTER_TABLE) \
		taxa.ifn=$(SET_TAX_LEGEND) \
		reps.ifn=$(SC_SUMMARY_UNIQUE) \
		res.ifn=$(SET_TAXA_RESOLVE) \
		fdir=$(TAXA_FDIR)/6_combined

# plot resolution
plot_taxa_resolve_barplot:
	$(_R) R/taxa_resolve_plot.r plot.barplot \
		ifn=$(SET_TAXA_RESOLVE) \
		fdir=$(TAXA_FDIR)/7_resolve

# resolved taxa
plot_taxa_resolve_details:
	$(_R) R/taxa_resolve_plot.r plot.details.resolve \
		ifn=$(SET_TAXA_RESOLVE) \
		fdir=$(TAXA_FDIR)/8_resolve_details

# nearest species taxa
plot_taxa_ref_details:
	$(_R) R/taxa_resolve_plot.r plot.details.ref \
		ifn=$(SET_TAXA_REPS) \
		fdir=$(TAXA_FDIR)/9_ref_details

make_taxa_plots: \
	copy_taxa_table plot_taxa_uniref plot_taxa_path \
	plot_taxa_trees \
	plot_taxa_limited_trees \
	plot_taxa_legend \
	plot_taxa_rep_detail \
	plot_taxa_rep plot_taxa_rep_analysis  \
	plot_taxa_single_tree \
	plot_taxa_summary_basic plot_taxa_summary_taxa \
	plot_taxa_resolve_barplot plot_taxa_resolve_details plot_taxa_ref_details

###########################################################
# obsolete
###########################################################

# plot all reps, with divergence and gene addition
plot_taxa_summary_old:
	$(_R) R/plot_taxa_reps.r plot.rep.summary \
		ifn.rep=$(SC_SUMMARY_UNIQUE) \
		ifn.legend=$(SET_TAXA_REP_LEGEND) \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		fdir=$(TAXA_FDIR)/summary

# shared genes have weaker matches in uniprot database
plot_taxa_shared_analysis:
	$(_R) R/taxa_shared_analysis.r plot.shared.analysis  \
		ifn.uniref=$(UNIREF_GENE_TAX_TABLE) \
		ifn.genes=$(TAXA_SET_GENES) \
		fdir=$(TAXA_FDIR)/shared_analysis

