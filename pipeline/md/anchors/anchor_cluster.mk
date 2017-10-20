
$(ANCHOR_CLUSTER_TABLE):
	@$(MAKE) make_set_clusters \
		SET_CLUSTER_DIR=$(ANCHOR_CLUSTER_DIR) \
		SET_CLUSTER_INPUT_MATRIX=$(ANCHOR_CLUSTER_INPUT_MATRIX) \
		SET_CLUSTER_INPUT_FIELD=$(ANCHOR_CLUSTER_INPUT_FIELD) \
		SET_CLUSTER_INPUT_MAP=$(ANCHOR_CLUSTER_INPUT_MAP) \
		SET_CLUSTER_TABLE=$(ANCHOR_CLUSTER_TABLE_LOCAL) \
		SET_CLUSTER_PREFIX=a \
		SET_GENES_TABLE=$(GENE_TABLE) \
		SET_CLUSTER_QSUB_DIR=$(ANCHOR_QSUB_MAP_DIR) \
		SET_SELF=T \
		SET_TITLE=predicted
	@mkdir -p $(ORDER_DIR)
	cp $(ANCHOR_CLUSTER_TABLE_LOCAL) $(ANCHOR_CLUSTER_TABLE)
make_anchor_cluster: $(ANCHOR_CLUSTER_TABLE)

plot_anchor_cluster: $(ANCHOR_CLUSTER_TABLE)
	@$(MAKE) plot_set_clusters \
		SET_CLUSTER_DIR=$(ANCHOR_CLUSTER_DIR) \
		SET_CLUSTER_INPUT_MATRIX=$(ANCHOR_CLUSTER_INPUT_MATRIX) \
		SET_CLUSTER_INPUT_FIELD=$(ANCHOR_CLUSTER_INPUT_FIELD) \
		SET_CLUSTER_INPUT_MAP=$(ANCHOR_CLUSTER_INPUT_MAP) \
		SET_CLUSTER_TABLE=$(ANCHOR_CLUSTER_TABLE) \
		SET_TITLE=predicted \
		SET_GENES_TABLE=$(GENE_TABLE) \
		SET_CLUSTER_QSUB_DIR=$(ANCHOR_QSUB_MAP_DIR) \
		SET_SELF=T \
		SET_CLUSTER_FDIR=$(ANCHOR_CLUSTER_FDIR)

.PHONY: make_anchor_cluster plot_anchor_cluster
