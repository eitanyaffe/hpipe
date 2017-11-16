SIMULATED_ID?=EC_BT_mix

# assembly
ASSEMBLY_COORDS?=$(ADIR)/coords
ASSEMBLY_READS1?=$(ASSEMBLY_READ_DIR)/R1.fastq
ASSEMBLY_READS2?=$(ASSEMBLY_READ_DIR)/R2.fastq
ASSEMBLY_READS_DONE?=$(ADIR)/.done

# hi-c
HIC_COORDS?=$(HIC_DIR)/hic_coords
HIC_READS1?=$(HIC_READ_DIR)/R1.fastq
HIC_READS2?=$(HIC_READ_DIR)/R2.fastq
HIC_READS_DONE?=$(HIC_DIR)/.done

##################################################################################################
# assembly
##################################################################################################

$(ASSEMBLY_COORDS): $(REF_GENOME_TABLE)
	$(call _start,$(ADIR))
	$(call _time,$(ADIR)) $(_R) $(_md)/R/simulate.r generate.assembly.read.coords.table \
		genome.table.ifn=$(REF_GENOME_TABLE) \
		xcoverage.base=$(BASE_XCOVERAGE) \
		read.length=$(READ_LENGTH) \
		insert=$(INSERT) \
		insert.sd=$(INSERT_SD) \
		base.length=$(BASE_MOLECULE_LENGTH) \
		ofn=$@
	$(_end)
assembly_coords: $(ASSEMBLY_COORDS)

$(ASSEMBLY_READS_DONE): $(ASSEMBLY_COORDS)
	$(call _start,$(ASSEMBLY_READ_DIR))
	$(call _time,$(ASSEMBLY_READ_DIR)) perl $(_md)/pl/generate_reads_pair.pl \
		$(ASSEMBLY_COORDS) \
		$(REF_FASTA) \
		$(ASSEMBLY_READS1) \
		$(ASSEMBLY_READS2)
	$(_end_touch)
assembly_reads: $(ASSEMBLY_READS_DONE)

##################################################################################################
# Hi-C
##################################################################################################

$(HIC_COORDS): $(REF_GENOME_TABLE)
	$(call _start,$(HIC_DIR))
	$(call _time,$(HIC_DIR)) $(_R) $(_md)/R/simulate.r generate.hic.read.coords.table \
		genome.table.ifn=$(REF_GENOME_TABLE) \
		nreads=$(HIC_READ_COUNT) \
		read.length=$(READ_LENGTH) \
		exponent=$(HIC_POWER_EXPONENT) \
		p.noise=$(HIC_NOISE_PROB) \
		p.distal=$(HIC_DISTAL_PROB) \
		ofn=$@
	$(_end)
hic_coords: $(HIC_COORDS)

$(HIC_READS_DONE): $(HIC_COORDS)
	$(call _start,$(HIC_READ_DIR))
	$(call _time,$(HIC_READ_DIR)) perl $(_md)/pl/generate_reads_pair.pl \
		$(HIC_COORDS) \
		$(REF_FASTA) \
		$(HIC_READS1) \
		$(HIC_READS2)
	$(_end_touch)
hic_reads: $(HIC_READS_DONE)

##################################################################################################
# Main
##################################################################################################

sim: $(ASSEMBLY_READS_DONE) $(HIC_READS_DONE)
