####################################################################
# Sup Figures
####################################################################

plot_timeline:
	$(_R) R/plot_timeline.r plot.timeline \
		ids=$(TDATASETS) \
		disturb.days=$(TDISTURB_DAYS) \
		labels=$(TLABELS) \
		days=$(TDAYS) \
		zoom.days=$(TZOOM_DAYS) \
		fdir=$(RES_FDIR)/timeline

plot_hic_vs_sg_abundance:
	$(_R) R/hic_vs_sg_abundance.r plot.scatter \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.anchors=$(ANCHOR_CLUSTER_TABLE) \
		ifn.hic.coverage=$(HIC_COVERAGE_TABLE) \
		ifn.sg.coverage=$(SG_COVERAGE_TABLE) \
		ifn.taxa.legend=$(SET_TAX_LEGEND) \
		fdir=$(RES_FDIR)/hic_vs_sg_abundance
print_response:
	$(_R) R/response_info.r print.info \
		ifn.info=$(RESPONSE_INFO) \
		odir=$(RES_FDIR)

plot_response_anchor_distrib:
	$(_R) R/response_plot.r plot.response.anchor.distrib \
		ifn.obs=$(RESPONSE_CONTIG_OBSERVED) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		fdir=$(RES_FDIR)/hierarchy_leaf_height_distrib

# response correlation matrix per anchors/clouds
plot_response_matrices:
	$(_R) R/response_anchor.r plot.response.matrix.anchor \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.norm=$(RESPONSE_CONTIG_NORM) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		fdir=$(RES_FDIR)/linkage

plot_anchor_control:
	$(_R) R/response_plot.r plot.anchor.control \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.contigs=$(CONTIG_TABLE) \
		ifn.mean=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.norm=$(RESPONSE_CONTIG_NORM) \
		fdir=$(RES_FDIR)/anchor_control

response_plot_base: plot_timeline print_response plot_response_matrices plot_anchor_control plot_response_anchor_distrib

####################################################################
# Host response
####################################################################

plot_sample_matrix:
	$(_R) R/sample_matrix_plot.r plot.sample.matrix \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.response=$(RESPONSE_CONTIG_NORM) \
		ifn.contigs=$(ANCHOR_TABLE) \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		labels=$(TLABELS) \
		fdir=$(RES_FDIR)/sample_matrix

plot_response_anchors:
	$(_R) R/response_plot.r plot.anchor.patterns \
		ifn.order=$(ANCHOR_CLUSTER_TABLE) \
		ifn.majors=$(RESPONSE_MAJORS) \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.top95=$(RESPONSE_PATTERN_TOP95) \
		ifn.top75=$(RESPONSE_PATTERN_TOP75) \
		ifn.bottom05=$(RESPONSE_PATTERN_BOTTOM05) \
		ifn.bottom25=$(RESPONSE_PATTERN_BOTTOM25) \
		ifn.taxa=$(SET_TAXA_REPS) \
		ifn.taxa.legend=$(SET_TAX_LEGEND) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		base.ids=$(RESPONSE_BASE_IDS) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		labels=$(TLABELS) \
		fdir=$(RES_FDIR)/host_response

plot_response_genus:
	$(_R) R/response_plot.r plot.genus.patterns \
		ifn.rep=$(SET_TAXA_REPS) \
		ifn.majors=$(RESPONSE_MAJORS) \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.taxa.legend=$(SET_TAX_LEGEND) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		labels=$(TLABELS) \
		ifn.taxa=$(SET_TAXA_TABLE) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		fdir=$(RES_FDIR)/genus_response

plot_response_group:
	$(_R) R/response_plot.r plot.group.patterns \
		ifn.rep=$(SET_TAXA_REPS) \
		ifn.majors=$(RESPONSE_MAJORS) \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.detection=$(RESPONSE_CONTIG_NORM_DETECTION) \
		ifn.groups=$(GROUP_ANCHOR_TABLE) \
		labels=$(TLABELS) \
		disturb.ids=$(RESPONSE_DISTURB_IDS) \
		fdir=$(RES_FDIR)/group_response

plot_response_classes:
	$(_R) R/response_plot.r plot.classes \
		ifn.median=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.class=$(RESPONSE_ANCHOR_ORDER) \
		fdir=$(RES_FDIR)/class_response

plot_response_vs_identity:
	$(_R) R/response_plot.r plot.response.vs.identity \
		ifn.majors=$(RESPONSE_MAJORS) \
		ifn.mean=$(RESPONSE_PATTERN_MEDIAN) \
		ifn.identity=$(ANCHOR_MATRIX_TABLE) \
		base.ids=$(RESPONSE_BASE_IDS) \
		fdir=$(RES_FDIR)/identity

response_plot_host: plot_sample_matrix plot_response_anchors plot_response_vs_identity plot_response_genus plot_response_group

####################################################################
# Experimental plots
####################################################################

# scatter of CA contact enrichment vs CA temporal correlation
plot_contact_vs_temporal:
	$(_R) R/response_plot.r plot.contact.vs.temporal \
		ifn.contigs=$(CONTIG_TABLE) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.matrix=$(CA_MATRIX) \
		min.contacts=$(CA_MIN_CONTACTS) \
		ifn.norm=$(RESPONSE_CONTIG_NORM) \
		ifn.anchors=$(RESPONSE_ANCHOR_GLOBAL) \
		fdir=$(RES_FDIR)/contact_vs_temporal

####################################################################
# obsolete plots
####################################################################

# complete matrix, colored by anchors. Cannot be plotted since we cluster without the anchors
obs_plot_response_raw:
	$(_R) R/response_plot.r plot.response.matrix \
		ifn.obs=$(RESPONSE_CONTIG_OBSERVED) \
		ifn.ca=$(CA_ANCHOR_CONTIGS) \
		ifn.cluster=$(RESPONSE_CONTIG_CLUSTER) \
		ifn.order=$(RESPONSE_CONTIG_ORDER) \
		fdir=$(RES_FDIR)/raw_cluster_analysis
