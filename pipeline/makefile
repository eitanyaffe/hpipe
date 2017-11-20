# makeshift
include makeshift/makeshift.mk

#####################################################################################################
# macros
#####################################################################################################

# count lines matching pattern path/*R1* and path/*R2*
# 1: append to this output file two lines: (field_R1,count1) and (field_R2,count2)
# 2: field
# 3: path
# 4: lines per read
_lc=echo -e $2_R1'\t'`cat $3/*R1* | wc -l`'\t'$4 >> $1; echo -e $2_R2'\t'`cat $3/*R2* | wc -l`'\t'$4 >> $1
_lc_pair=echo -e $2'\t'`cat $3/P* | wc -l`'\t'$4 >> $1

#####################################################################################################
# classes
#####################################################################################################

# for preproc
$(call _class,library,LIB_ID LIB_INPUT_DIRS PREPROC_MODES)

# for making assemblies
$(call _class,assembly,ASSEMBLY_ID ASSEMBLY_LIB_IDS)

# mapping fastq to assemblies
$(call _class,map,ASSEMBLY_ID LIB_ID)
$(call _class,map2,ASSEMBLY_ID LIB_ID MAP_INPUT)

# single end mapping
$(call _class,se_map,ASSEMBLY_ID LIB_ID MAP_INPUT)

# anchors
$(call _class,anchor,ANCHOR ASSEMBLY_ID LIB_ID)

# assembly+anchor+library
$(call _class,ca,ASSEMBLY_ID ANCHOR LIB_ID)

# assembly+anchor+library
$(call _class,response,ASSEMBLY_ID ANCHOR LIB_ID RESPONSE_LIB_ID)

# loop over classes
map_loop:
	@$(MAKE) class_loop class=map
se_map_loop:
	@$(MAKE) class_loop class=se_map
assembly_loop:
	@$(MAKE) class_loop class=assembly
lib_loop:
	@$(MAKE) class_loop class=library
ca_loop:
	@$(MAKE) class_loop class=ca
response_loop:
	@$(MAKE) class_loop class=response
anchor_loop:
	@$(MAKE) class_loop class=anchor


#####################################################################################################
# modules
#####################################################################################################

# add modules
$(call _module,global.mk)

# download and simulate maps
$(call _module,md/ref/ref_int.mk)

# library pre-process
$(call _module,md/preproc/preproc_int.mk)

# genome assembly
$(call _module,md/assembly/assembly_int.mk)

# mapping using bwa
$(call _module,md/map/map_int.mk)

# initial contig clustering
$(call _module,md/ccluster/ccluster_int.mk)

# anchors
$(call _module,md/anchors/anchors_int.mk)

# docker
$(call _module,md/docker/docker_int.mk)

# genes and anchor/genes
$(call _module,md/genes/genes_int.mk)

# uniref
$(call _module,md/taxa/taxa_int.mk)

# cell response
$(call _module,md/response/response_int.mk)

# export tables
$(call _module,md/export/export_int.mk)

#####################################################################################################
# config file
#####################################################################################################

# c: config file
# NOTE: customize after modules loaded since we use module variables for class functions
c?=config/ref/n5.cfg
include $(c)

$(call _set_user_title,config: $(c))

#####################################################################################################
# module / step
#####################################################################################################

# select module
m?=response
$(call _active_module,$(m))

# select step
$(call _step,$(s))

# $(info c=$(c), m=$(m), s=$(s))

#####################################################################################################
# compile all binaries
#####################################################################################################

init:
	@$(MAKE) m=preproc preproc_init
	@$(MAKE) m=ccluster ccluster_init
	@$(MAKE) m=anchors trim_init
	@$(MAKE) m=anchors model_init_bin
	@$(MAKE) m=response response_init
	@$(MAKE) m=taxa init_sc

#####################################################################################################
# pre-process fastq libraries
#####################################################################################################

# just library pre-process
lib:
	@$(MAKE) m=preproc make_preproc

# for library
plot_lib:
	@$(MAKE) m=preproc make_plot_preproc

#####################################################################################################
# make assembly
#####################################################################################################

# LIB_ID and PREPROC_MODE must be set
assembly:
	@$(MAKE) m=assembly make_assembly

# for assembly
plot_assembly:
	@$(MAKE) m=assembly make_plot_assembly

#####################################################################################################
# Map reads to assembly
#####################################################################################################

# map w/o pairing and compute coverage
map_se:
	@$(MAKE) m=map coverage

# map paired-end reads to assembly, and compute coverage
map_pe:
	@$(MAKE) m=map map_all

#	$(foreach LIB_ID,$(LIB_IDS),$(MAKE) map_pe LIB_ID=$(LIB_ID); $(ASSERT);)
# map libs
#map_assembly_libs:
#	@echo looping over maps: $(call _class_get_instances,map)
#	@$(MAKE) class_loop class=map_libs t=map_libs

# for map
plot_map:
	@$(MAKE) m=map make_plot_map_summary

#####################################################################################################
# Predict genes and blast over uniref
#####################################################################################################

# over assembly
genes:
	@$(MAKE) m=genes make_genes

uniref:
	@$(MAKE) m=genes make_uniref

#####################################################################################################
# contig matrix
#####################################################################################################

# filter out contacts by contig similarity and generate matrix
cmatrix:
	@$(MAKE) m=ccluster contig_matrix

cis_decay:
	@$(MAKE) m=ccluster make_cis_decay

# for assembly
plot_contig_matrix:
	@$(MAKE) m=ccluster compare_cc_matrix

# for map
plot_cis_decay:
	@$(MAKE) m=ccluster make_plot_cis_decay

# for assembly
plot_cis_decay_assembly:
	@$(MAKE) m=ccluster make_plot_cis_decay_assembly

#####################################################################################################
# Compute anchor set
#####################################################################################################

# anchor seeds
ianchors:
	@$(MAKE) m=ccluster make_ianchors

# models
model_init:
	@$(MAKE) m=anchors model MDL_STAGE=initial
model_final:
	@$(MAKE) m=anchors model MDL_STAGE=final
model_final_cellular:
	@$(MAKE) m=anchors model MDL_STAGE=final_cellular

# the whole thing
anchors:
	@$(MAKE) m=ccluster make_ianchors
	@$(MAKE) m=anchors coords2fends
	@$(MAKE) m=anchors model MDL_STAGE=initial
	@$(MAKE) m=anchors trim_anchors
	@$(MAKE) m=anchors model MDL_STAGE=final

# use existing anchor
use_anchors:
	@$(MAKE) m=anchors coords2fends
	@$(MAKE) m=anchors model MDL_STAGE=final

#####################################################################################################
# Compute gene/anchor matrix (ca)
#####################################################################################################

ca_matrix:
	@$(MAKE) m=anchors coords2fends
	@$(MAKE) m=anchors model MDL_STAGE=final
	@$(MAKE) m=anchors make_ca

#####################################################################################################
# Cluster anchors
#####################################################################################################

cluster:
	@$(MAKE) m=anchors make_anchor_cluster

#####################################################################################################
# Reference genome
#####################################################################################################

# over assembly
ref_genes:
	@$(MAKE) m=genes make_ref

# over anchor
ref_cluster:
	@$(MAKE) m=anchors make_ref_cluster make_ref_lookup

#####################################################################################################
# Compute anchor NCBI taxa trees
#####################################################################################################

# over anchors
taxa:
	@$(MAKE) m=taxa make_taxa

#####################################################################################################
# Temporal response
#####################################################################################################

# over response
response:
	@$(MAKE) m=response response_network

# over assembly
response_compare:
	@$(MAKE) m=response response_network_compare

#####################################################################################################
# Compute assembly props as a function of read depth
#####################################################################################################

rarify_reads:
	@$(MAKE) m=assembly plot_rarify

#####################################################################################################
#####################################################################################################

#####################################################################################################
# package
#####################################################################################################

BCC_VERSION=.version
BCC_FILES=.version makefile md global.mk input
pack:
	@$(MAKE) package PPATHS="$(BCC_FILES)" PVERSION=$(BCC_VERSION) PTITLE=bcc POUTDIR=versions

#####################################################################################################
# pipe combos
#####################################################################################################

# only libs
pp_libs:
	@$(MAKE) lib_loop t=lib

# basic pipe (assembly+mapping+cis_decay+contig_matrix)
pp_assembly:
	@$(MAKE) lib_loop t=lib
	@$(MAKE) assembly_loop t=assembly
	@$(MAKE) map_loop t=map_pe
	@$(MAKE) map_loop t=cis_decay
	@$(MAKE) map_loop t=cmatrix

# genes
pp_genes:
	@$(MAKE) assembly_loop t=genes

# blast to uniref
pp_uniref:
	@$(MAKE) assembly_loop t=uniref

# anchor inference
pp_anchors:
	@$(MAKE) anchor_loop t=anchors
	@$(MAKE) ca_loop t=ca_matrix
	@$(MAKE) anchor_loop t=cluster
	@$(MAKE) m=anchors anchor_loop t=make_anchor_groups
	@$(MAKE) anchor_loop m=anchors t=anchor_stats

# flat file structure
pp_output_basic:
	@$(MAKE) m=export anchor_loop t=flat_basic

pp_checkm:
	@$(MAKE) anchor_loop m=anchors t=make_checkm

pp_taxa:
	@$(MAKE) anchor_loop t=taxa

# single ended (no preproc) mapping, for temporal shot gun2
pp_map_se:
	@$(MAKE) se_map_loop t=map_se

# response
pp_response:
	@$(MAKE) response_loop m=response t=response
pp_response_compare:
	@$(MAKE) anchor_loop m=response t=response_compare

pp_all: pp_assembly pp_genes pp_uniref pp_anchors pp_checkm pp_taxa pp_map_se pp_response

# assembly, genes, anchor/union pairs
pp_basic: pp_assembly pp_genes pp_anchors pp_output_basic

# only up to taxa (no disturbance)
pp_all_short: pp_assembly pp_genes pp_anchors pp_checkm pp_taxa

# use existing anchor set
pp_use_anchors:
	@$(MAKE) lib_loop t=lib
	@$(MAKE) map_loop t=map_pe
	@$(MAKE) map_loop t=cis_decay
	@$(MAKE) map_loop t=cmatrix
	@$(MAKE) anchor_loop t=use_anchors
	@$(MAKE) ca_loop t=ca_matrix

#####################################################################################################
# plot combos
#####################################################################################################

# basic plots
pl_qc:
	@$(MAKE) lib_loop t=plot_lib
	@$(MAKE) assembly_loop t=plot_assembly
	@$(MAKE) assembly_loop t=plot_map
	@$(MAKE) map_loop t=plot_cis_decay
	@$(MAKE) assembly_loop t=plot_cis_decay_assembly
	@$(MAKE) assembly_loop t=plot_contig_matrix

pl_ccluster:
	@$(MAKE) m=ccluster anchor_loop t=iplots

# anchor related
pl_anchors:
	@$(MAKE) m=anchors anchor_loop t=contact_stats
	@$(MAKE) m=anchors anchor_loop t=plot_trim
	@$(MAKE) m=anchors anchor_loop t=plot_model MDL_STAGE=initial
	@$(MAKE) m=anchors anchor_loop t=plot_model MDL_STAGE=final
	@$(MAKE) m=anchors anchor_loop t=plot_model MDL_STAGE=final_cellular
	@$(MAKE) m=anchors anchor_loop t=plot_anchor_cluster
	@$(MAKE) m=anchors anchor_loop t=make_info_plots
	@$(MAKE) m=anchors ca_loop t=make_ca_plots
	@$(MAKE) m=anchors anchor_loop t=make_group_plots

# taxa related
pl_taxa:
	@$(MAKE) anchor_loop m=taxa t=make_taxa_plots

# checkm plots
pl_checkm:
	@$(MAKE) anchor_loop m=anchors t=make_checkm_plots

# cipro response
pl_response:
	@$(MAKE) anchor_loop m=response t=response_plot_base
	@$(MAKE) anchor_loop m=response t=response_plot_host
	@$(MAKE) response_loop m=response t=response_plot_network

pl_response_compare:
	@$(MAKE) anchor_loop m=response t=plot_network_compare

pl_all: pl_qc pl_ccluster pl_anchors pl_checkm pl_taxa pl_response

pl_basic: pl_qc pl_ccluster pl_anchors

# no response
pl_all_short: pl_qc pl_ccluster pl_anchors pl_checkm pl_taxa

#####################################################################################################
# Hi-C simulated data
#####################################################################################################

# simulate shotgun and Hi-C reads
pp_simulate:
	@$(MAKE) m=ref download
	@$(MAKE) m=ref sim

# extra for reference tests
pp_ref:
	@$(MAKE) assembly_loop t=ref_genes
	@$(MAKE) anchor_loop t=ref_cluster

pl_ref:
	@$(MAKE) anchor_loop m=genes t=plot_classify
	@$(MAKE) anchor_loop m=anchors t=plot_ref_cluster
	@$(MAKE) anchor_loop m=anchors t=make_ref_plots

#####################################################################################################
# Gene sharing patterns
#####################################################################################################

# compute gene sharing between isolate genomes
pp_gene_sharing:
	@$(MAKE) m=genes gr_all
	@$(MAKE) m=anchors gs_all

pp_docker_image:
	@$(MAKE) m=docker docker_build
