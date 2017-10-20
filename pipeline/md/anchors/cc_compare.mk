DATASET1=pre_lib_hic_simple
DATASET2=post_lib_hic_simple

plot_cc_compare:
	$(_R) R/cc_compare.r plot.cc.compare \
		ifn1=$(call reval,CA_MATRIX,DATASET=$(DATASET1)) \
		ifn2=$(call reval,CA_MATRIX,DATASET=$(DATASET2)) \
		fdir=$(SET_FIGURE_DIR)/cc_compare
