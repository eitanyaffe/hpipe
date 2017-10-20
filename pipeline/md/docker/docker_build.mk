
# manualy download MetaGeneMark from: http://topaz.gatech.edu/Genemark/license_download.cgi

RSYNC_CMD=rsync -rq
docker_files:
#	mkdir -p $(DOCKER_DIR)
	cp -p $(SRC_DOCKER_FILE) $(DOCKER_FILE)
#	cp -p input/download/gm_key_64.gz input/download/MetaGeneMark_linux_64.tar.gz $(DOCKER_DIR)
	touch $(DOCKER_DIR)/.is_container
#	$(RSYNC_CMD) md $(DOCKER_DIR)
#	$(RSYNC_CMD) input $(DOCKER_DIR)
#	$(RSYNC_CMD) -L makeshift $(DOCKER_DIR)
#	cp makefile $(DOCKER_DIR)/
#	cp global.mk $(DOCKER_DIR)/

docker_build: docker_files
	docker build -f $(DOCKER_FILE) -t $(DOCKER_ID) $(DOCKER_DIR)

docker_push:
	docker tag $(DOCKER_ID) $(DOCKER_NAME)/$(DOCKER_REPOSITORY):$(DOCKER_TAG)
	docker push $(DOCKER_NAME)/$(DOCKER_REPOSITORY):$(DOCKER_TAG)
