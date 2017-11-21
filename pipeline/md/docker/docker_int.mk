#####################################################################################################
# register module
#####################################################################################################

units=docker_build.mk
$(call _register_module,docker,$(units),docker,)

#####################################################################################################
# docker
#####################################################################################################

DOCKER_ID?=hpipe
#DOCKER_ID?=centos_r

BASE_DOCKER_DIR?=$(BASE_OUTDIR)/docker_images
DOCKER_DIR?=$(BASE_DOCKER_DIR)/$(DOCKER_ID)
SRC_DOCKER_FILE?=input/docker/$(DOCKER_ID)
DOCKER_FILE?=$(DOCKER_DIR)/Dockerfile

DOCKER_NAME=eitanyaffe
DOCKER_REPOSITORY=hpipe
DOCKER_TAG=v1.01
