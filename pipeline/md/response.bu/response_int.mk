#####################################################################################################
# register module
#####################################################################################################

pr_vars:=
units:=response.mk
$(call _register_module,response,$(units),global,$(pr_vars))

#####################################################################################################

RESPONSE_ID?=cipro

RESPONSE_DIR?=$(ASSEMBLY_DIR)/response/$(RESPONSE_ID)
RESPONSE_CONTIG_OBSERVED?=$(RESPONSE_DIR)/contigs_observed
RESPONSE_CONTIG_EXPECTED?=$(RESPONSE_DIR)/contigs_expected
RESPONSE_CONTIG_NORM?=$(RESPONSE_DIR)/contigs_norm

RESPONSE_TABLE?=$(RESPONSE_DIR)/table
RESPONSE_ORDER?=$(RESPONSE_DIR)/order

########################################################################
# cutoff clustering
########################################################################

RESPONSE_CONTIG_TABLE_NORM?=$(RESPONSE_DIR)/contigs_n
RESPONSE_CLUSTER_THRESHOLD?=0.98
RESPONSE_CONTIG_CLUSTER_METHOD?=pearson
RESPONSE_CLUSTER_ID?=$(RESPONSE_CLUSTER_THRESHOLD)_$(RESPONSE_CONTIG_CLUSTER_METHOD)
RESPONSE_CONTIG_CLUSTERS?=$(RESPONSE_DIR)/contigs_clustered_$(RESPONSE_CLUSTER_ID)
RESPONSE_CLUSTER_HIST=$(RESPONSE_DIR)/hist_$(RESPONSE_CLUSTER_ID)

########################################################################
# kmeans clustering
########################################################################

RESPONSE_CONTIG_TABLE_FOLD?=$(RESPONSE_DIR)/contigs_fold
RESPONSE_KMEANS_MAX_FOLD?=100
RESPONSE_KMEANS_MIN_CONTIG_COVERAGE?=30

# random or diameter
RESPONSE_KMEANS_SEED_METHOD?=random

RESPONSE_KMEANS_N_CLUSTERS?=100
RESPONSE_KMEANS_SEED?=1010
RESPONSE_KMEANS_PREFIX?=$(RESPONSE_DIR)/kmeans
RESPONSE_KMEANS_CONTIGS?=$(RESPONSE_KMEANS_PREFIX).clusters
RESPONSE_KMEANS_CENTROID_MEAN?=$(RESPONSE_KMEANS_PREFIX).centroid_mean
RESPONSE_KMEANS_CENTROID_SD?=$(RESPONSE_KMEANS_PREFIX).centroid_sd

# order centroids by hclust
RESPONSE_KMEANS_ORDER?=$(RESPONSE_KMEANS_PREFIX).order

########################################################################

TDATASETS=cipro_S1 cipro_S2 cipro_S3 cipro_S4 cipro_S5 cipro_S6 cipro_S7 cipro_S8 cipro_S9 cipro_S10 cipro_S11 cipro_S12

