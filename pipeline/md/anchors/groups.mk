
#####################################################################################################
# partition anchors into groups
#####################################################################################################

GROUP_ANCHORS_DONE?=$(GROUP_DIR)/.done_groups
$(GROUP_ANCHORS_DONE):
	$(call _start,$(GROUP_DIR))
	$(_R) R/anchor_groups.r group.anchors \
		mat.ifn=$(ANCHOR_MATRIX_TABLE) \
		min.identity=$(MEAN_IDENTITY_THRESHOLD) \
		group.threshold=$(GROUP_THRESHOLD) \
		order.ifn=$(ANCHOR_CLUSTER_TABLE) \
		group.ofn=$(GROUP_TABLE) \
		anchor.ofn=$(GROUP_ANCHOR_TABLE)
	$(_end_touch)

make_anchor_groups: $(GROUP_ANCHORS_DONE)

#####################################################################################################
# plotting
#####################################################################################################

plot_group_dendrogram:
	$(_start)
	$(_R) R/anchor_groups.r plot.group.dendrogram \
		mat.ifn=$(ANCHOR_MATRIX_TABLE) \
		min.identity=$(MEAN_IDENTITY_THRESHOLD) \
		group.threshold=$(GROUP_THRESHOLD) \
		order.ifn=$(ANCHOR_CLUSTER_TABLE) \
		fdir=$(GROUP_FDIR)
	$(_end)

make_group_plots: plot_group_dendrogram
