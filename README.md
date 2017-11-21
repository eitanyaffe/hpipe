# hpipe: A tool for the analysis of metagenomic Hi-C data

The hpipe tool processes metagenomic Hi-C maps to infer a set of underlying
genome anchor/union pairs. Each *genome anchor* is a collection contigs that
are part of the genomes of one or more existing strains, while each genome
anchor has a matching *genome union*, that is the combined genomes of
those strains.

As a service, the pipeline starts of from raw reads and performs various
initial steps including quality trimming, removal of adapter sequences
and metagenomic assembly. The core functionality of the pipeline is
(1) to infer a background model of spurious contacts between contigs that
belong to different cells, and (2) to infer a set of anchors.

The tool was developed by Eitan Yaffe, at David Relman's lab, Stanford
University School of Medicine. It is distributed under the GNU General
Public License v3.0. If you have questions or comments please contact Eitan
Yaffe at eitan.yaffe@gmail.com.

## Pre-requisites

The tool requires docker (https://www.docker.com), perl and bash. All work is
performed within a docker container, simplifying the installation and making the
tool compatible with most linux systems. The current version was tested
under CentoOS 6.9.

One optional pre-requisite is DeconSeq (http://deconseq.sourceforge.net/),
a pipeline for the removal of human reads. By default the tool does not remove
human reads.

## Installation

1. Select a working directory. Here we use /work/hpipe as an example.
```
mkdir -p /work
cd /work
```

2. Get source code from github.
```
git clone https://github.com/eitanyaffe/hpipe.git
```

3. Set the environment variable HPIPE_DIR to point to the hpipe dir.
For example, if using bash, you can add the following line to your .bashrc:
```
export HPIPE_DIR=/work/hpipe
```

4. Add the script `$HPIPE_DIR/hpipe.pl` to your path. To do that either copy the
script to some common directory (e.g. /usr/local/bin), or add $HPIPE_DIR to your
path.

## Quick Start and Usage

To verify hpipe has been successfully installed run the following commands:

1. Start an hpipe container. Each project requires a separate container.
```
./hpipe.pl start -c config/template/basic.cfg
```

2. Infer anchor-union pairs. This should take less than 10 minutes for the test dataset.
```
./hpipe.pl run -c config/template/basic.cfg
```

3. When done close the container.
```
./hpipe.pl stop -c config/template/basic.cfg
```

Output files should be generated under projects/project1/output.

To run hpipe create a configuration directory based of the original template or
some other pre-existing template. For example:
```
cp -r config/template config/my_project
```

Then edit the various files inside the configuration directory (see details below)
and run as above.

## Input and Output

Input and output files and directories are defined using a single configuration
directory per each project.  You can see an example of a configuration directory
under config/template. To start a new project please copy config/template
and modify the various files there.

The configuration directory contains 3 files:
* project_id: Unique string identifier of project (`PROJECT_ID`)
* path_vars: Paths to input and output directories.
* user_vars: Customizable parameters.
* basic.cfg: Experimental design. (Typically, do not edit this file.)

### Input

The pipeline uses as input two raw reads of DNA libraries.

* Assembly library, used to generate a metagenomic assembly. The location of the
fastq file is determined by the `ASSEMBLY_INPUTDIR` in the path_vars file.

* Hi-C library, used to infer the background model and infer anchor/genome
pairs. The location of the fastq file is determined by the `HIC_INPUTDIR` in
the path_vars file.

User-defined parameters in user_vars include:
* CUTTER_TITLE: The restriction enzyme used.
* CUTTER_SITE: The cutter site.
* REMOVE_HUMAN: [T|F] Should remove human reads using Deconseq

Each project is associated with a docker container. The container is setup
to point to the various paths exposed by the user in the path_vars file.

### Output

The pipeline outputs all files under `BASE_OUTDIR`, defined in path_vars.
The output of each project is placed under `BASE_OUTDIR/PROJECT_ID/output`.
Temporary files are placed under `BASE_TMPDIR`.

Output files include:

* contig.fasta: fasta file of the metagenomic assembly files.
* contig.table: table with all contigs and their length (in bp)
* anchor.table: table with all genome anchors. Fields:
  * contig: contig identifier.
  * anchor: anchor identifier.
* contig_anchor.table: table with all genome unions. Fields:
  * contig: contig identifier.
  * anchor: anchor identifier.
  * observed: number of contacts between contig and anchor.
  * expected: number of expected spurious contacts between contig and anchor.
  * score: enrichment of observed the expected, in log_10.
* model/fend.table: table with fragment ends. Fields:
  * fend: fragment end identifier
  * frag: fragment identifier
  * strand: strand of fragment end
  * contig: contig identifier
  * coord: coordinate of fragment end
  * frag_len: length of fragment
  * abundance: abundance of contig
  * anchor: anchor indentifier, if none equals 0
  * abundance_bin: abundance bin
  * frag_len_bin: fragment length bin
* model/fragment_length.f, model/abundance.f: correction matrices
* model/fragment_length.bins, model/abundance.bins: ranges of the bins
* model/prior: constant prior probability value
