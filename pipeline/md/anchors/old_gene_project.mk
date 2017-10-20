
$(GENE_PROJECT_MAP):
	$(call _start,$(GENE_PROJECT_DIR))
	perl $(_md)/pl/gene_project_map.pl \
		$(REF_MAP_INPUT) \
		$(GENE_PROJECT_IDENTITY_THRESHOLD) \
		$(GENE_PROJECT_COVERAGE_THRESHOLD) \
		$@
	$(_end)

GENE_PROJECT_FWD_DONE?=$(GENE_PROJECT_DIR)/.done_fwd
$(GENE_PROJECT_FWD_DONE): $(GENE_PROJECT_MAP)
	$(_start)
	perl $(_md)/pl/gene_project.pl \
		$(REF_SET_INPUT1) \
		$(REF_FIELD_INPUT1) \
		$(REF_MAP_INPUT) \
		gene1 gene2 \
		$(GENE_PROJECT_IDENTITY_THRESHOLD) \
		$(GENE_PROJECT_COVERAGE_THRESHOLD) \
		$(GENE_PROJECT_FWD_SOURCE) \
		$(GENE_PROJECT_FWD_TARGET)
	$(_end_touch)

GENE_PROJECT_BCK_DONE?=$(GENE_PROJECT_DIR)/.done_bck
$(GENE_PROJECT_BCK_DONE): $(GENE_PROJECT_MAP)
	$(_start)
	perl $(_md)/pl/gene_project.pl \
		$(REF_SET_INPUT2) \
		$(REF_FIELD_INPUT2) \
		$(REF_MAP_INPUT) \
		gene2 gene1 \
		$(GENE_PROJECT_IDENTITY_THRESHOLD) \
		$(GENE_PROJECT_COVERAGE_THRESHOLD) \
		$(GENE_PROJECT_BCK_SOURCE) \
		$(GENE_PROJECT_BCK_TARGET)
	$(_end_touch)

plot_source_stats: $(GENE_PROJECT_FWD_DONE) $(GENE_PROJECT_BCK_DONE)
	$(_start)
	$(_R) R/gene_project.r plot.source.stats \
		order.ifn1=$(REF_CLUSTER_TABLE1) \
		order.ifn2=$(REF_CLUSTER_TABLE2) \
		fwd.source.ifn=$(GENE_PROJECT_FWD_SOURCE) \
		bck.source.ifn=$(GENE_PROJECT_BCK_SOURCE) \
		fdir=$(GENE_PROJECT_FDIR)
	$(_end)

####################################################################################
# summary
####################################################################################

SET_FWD_SUMMARY_TABLE_DONE?=$(GENE_PROJECT_DIR)/.done_fwd_summary
$(SET_FWD_SUMMARY_TABLE_DONE): $(GENE_PROJECT_FWD_DONE) $(GENE_PROJECT_BCK_DONE) $(SET_MATCH_DONE)
	$(_start)
	$(_R) R/gene_project.r project.sets \
		map.ifn=$(GENE_PROJECT_MAP) \
		source.field=gene1 target.field=gene2 \
		marked.ifn=$(GENE_PROJECT_FWD_TARGET) \
		set.ifn=$(REF_SET_INPUT2) \
		set.field=$(REF_FIELD_INPUT2) \
		set.other.ifn=$(REF_SET_INPUT1) \
		set.other.field=$(REF_FIELD_INPUT1) \
		match.ifn=$(SET_MATCH_TABLE1) \
		ofn.summary=$(SET_FWD_SUMMARY_TABLE) \
		ofn.report=$(SET_FWD_REPORT_TABLE) \
		fwd.flag=T
	$(_end_touch)

SET_BCK_SUMMARY_TABLE_DONE?=$(GENE_PROJECT_DIR)/.done_bck_summary
$(SET_BCK_SUMMARY_TABLE_DONE): $(GENE_PROJECT_FWD_DONE) $(GENE_PROJECT_BCK_DONE) $(SET_MATCH_DONE)
	$(_start)
	$(_R) R/gene_project.r project.sets \
		map.ifn=$(GENE_PROJECT_MAP) \
		source.field=gene2 target.field=gene1 \
		marked.ifn=$(GENE_PROJECT_BCK_TARGET) \
		set.ifn=$(REF_SET_INPUT1) \
		set.field=$(REF_FIELD_INPUT1) \
		set.other.ifn=$(REF_SET_INPUT2) \
		set.other.field=$(REF_FIELD_INPUT2) \
		match.ifn=$(SET_MATCH_TABLE2) \
		ofn.summary=$(SET_BCK_SUMMARY_TABLE) \
		ofn.report=$(SET_BCK_REPORT_TABLE) \
		fwd.flag=F
	$(_end_touch)

plot_target_stats: $(GENE_PROJECT_FWD_DONE) $(GENE_PROJECT_BCK_DONE) $(SET_FWD_SUMMARY_TABLE_DONE) $(SET_BCK_SUMMARY_TABLE_DONE)
	$(_start)
	$(_R) R/gene_project.r plot.target.stats \
		order.ifn1=$(REF_CLUSTER_TABLE1) \
		order.ifn2=$(REF_CLUSTER_TABLE2) \
		fwd.ifn=$(SET_FWD_SUMMARY_TABLE) \
		bck.ifn=$(SET_BCK_SUMMARY_TABLE) \
		fdir=$(GENE_PROJECT_FDIR)
	$(_end)

project_plots: plot_source_stats plot_target_stats
