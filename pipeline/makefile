# makeshift
include $(MAKESHIFT_ROOT)/makeshift-core/makeshift.mk

# for unit tests, run on locally installed modules
include tests.mk

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
$(call _class,assembly,ASSEMBLER ASSEMBLY_ID ASSEMBLY_LIB_IDS)

# mapping fastq to assemblies
$(call _class,map,ASSEMBLY_ID LIB_ID)
$(call _class,map2,ASSEMBLY_ID LIB_ID MAP_INPUT)

# anchors
$(call _class,anchor,ANCHOR ASSEMBLY_ID LIB_ID)

# assembly+anchor+library
$(call _class,ca,ASSEMBLY_ID ANCHOR LIB_ID)

# loop over classes
map_loop:
	@$(MAKE) class_loop class=map
assembly_loop:
	@$(MAKE) class_loop class=assembly
lib_loop:
	@$(MAKE) class_loop class=library
ca_loop:
	@$(MAKE) class_loop class=ca
anchor_loop:
	@$(MAKE) class_loop class=anchor

#####################################################################################################
# modules
#####################################################################################################

# add modules
$(call _module,global.mk)

# library pre-process
$(call _module_root,hpipe_preproc)

# genome assembly
$(call _module_root,hpipe_assembly)

# mapping using bwa
$(call _module_root,hpipe_map)

# initial contig clustering
$(call _module_root,hpipe_ccluster)

# anchors
$(call _module_root,hpipe_anchors)

# genes and anchor/genes
$(call _module_root,hpipe_genes)

# checkm
$(call _module_root,hpipe_checkm)

# uniref
$(call _module_root,hpipe_taxa)

# export tables
$(call _module_root,hpipe_export)

#####################################################################################################
# config file
#####################################################################################################

# c: config file
# NOTE: customize after modules loaded since we use module variables for class functions
c?=
include $(c)

$(call _set_user_title,config: $(c))

#####################################################################################################
# module / step
#####################################################################################################

# select module
m?=anchors
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

#####################################################################################################
# pre-process fastq libraries
#####################################################################################################

# just library pre-process
lib:
	@$(MAKE) m=preproc make_preproc

#####################################################################################################
# make assembly
#####################################################################################################

# LIB_ID and PREPROC_MODE must be set
assembly:
	@$(MAKE) m=assembly make_assembly

#####################################################################################################
# Map reads to assembly
#####################################################################################################

# map paired-end reads to assembly, and compute coverage
map_pe:
	@$(MAKE) m=map map_all

#####################################################################################################
# contig matrix
#####################################################################################################

# filter out contacts by contig similarity and generate matrix
cmatrix:
	@$(MAKE) m=ccluster contig_matrix

#####################################################################################################
# Compute anchor set
#####################################################################################################

anchors:
	@$(MAKE) m=ccluster make_ianchors
	@$(MAKE) m=anchors coords2fends
	@$(MAKE) m=anchors model MDL_STAGE=initial
	@$(MAKE) m=anchors trim_anchors
	@$(MAKE) m=anchors model MDL_STAGE=final
	@$(MAKE) m=anchors cell_mdl MDL_STAGE=final_cellular
	@$(MAKE) m=anchors model_cell MDL_STAGE=final_cellular

#####################################################################################################
# Compute contig/anchor matrix
#####################################################################################################

ca_matrix_basic:
	@$(MAKE) m=anchors coords2fends
	@$(MAKE) m=anchors model MDL_STAGE=final
	@$(MAKE) m=anchors make_ca_basic

#####################################################################################################
# pipe combos
#####################################################################################################

# basic pipe (assembly+mapping+cis_decay+contig_matrix)
pp_assembly:
	@$(MAKE) lib_loop t=lib
	@$(MAKE) assembly_loop t=assembly
	@$(MAKE) map_loop t=map_pe
	@$(MAKE) map_loop t=cmatrix

# anchor inference
pp_anchor_basic:
	@$(MAKE) anchor_loop t=anchors
	@$(MAKE) ca_loop t=ca_matrix_basic

# flat file structure
pp_output_basic:
	@$(MAKE) m=export anchor_loop t=hpipe_basic

# assembly and anchor/union pairs only
pp_basic: pp_assembly pp_anchor_basic pp_output_basic

