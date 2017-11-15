# utility scripts
PARSE_SHOWCOORD?=$(_md)/pl/show_coords_parse.pl
MUMMER_PARSE?=$(_md)/pl/mum_parse.pl

SIM_SPLIT_DONE?=$(SIM_SPLIT_DIR)/.done
$(SIM_SPLIT_DONE):
	$(call _start,$(SIM_SPLIT_DIR))
	$(_md)/pl/split_contigs.pl \
		$(CONTIG_TABLE) \
		$(CONTIG_FILE) \
		$(SIM_CHUNKS) \
		$(SIM_SPLIT_DIR)
	$(_end_touch)

$(SIM_DONE): $(SIM_SPLIT_DONE)
	$(call _start,$(SIM_RESULT_DIR))
	$(call _file_exists,$(MUMMER))
	$(_R) $(_md)/R/distrib_similarity.r distrib.similarity \
		mummer=$(MUMMER) \
		mummer.parse=$(MUMMER_PARSE) \
		split.dir=$(SIM_SPLIT_DIR) \
		odir=$(SIM_RESULT_DIR) \
		qsub.dir=$(SIM_QSUB_DIR) \
		batch.max.jobs=$(MAX_JOBS) \
		total.max.jobs.fn=$(MAX_JOBS_FN) \
		dtype=$(DTYPE) \
		show.coords=$(SHOWCOORD) \
		parse.show.coords.script=$(PARSE_SHOWCOORD) \
		jobname=contig_sim
	$(_end_touch)
csim: $(SIM_DONE)
