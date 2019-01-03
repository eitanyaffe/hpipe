# hpipe: A tool for the analysis of metagenomic Hi-C data

The hpipe tool processes metagenomic Hi-C maps to infer a set of underlying
genome anchor/union pairs. Each *genome anchor* is a collection contigs that
are part of the genomes of one or more residing strains, while each genome
anchor has a matching *genome union*, that is the combined genomes of
those strains.

The pipeline recieves as input a shotgun libary and a Hi-C library,
and performs various initial steps including quality trimming, removal of
adapter sequences and metagenomic assembly. The core functionality of the
pipeline is (1) to infer a background model of spurious contacts between
contigs that belong to different cells, and (2) to infer a set of anchors.

The tool was developed by Eitan Yaffe, at David Relman's lab, Stanford
University School of Medicine. It is distributed under the GNU General
Public License v3.0. If you have questions or comments please contact Eitan
Yaffe at eitan.yaffe@gmail.com.

## Pre-requisites

The tool requires docker, perl and bash. All work is performed within a
docker container, simplifying the installation and making the tool
compatible with most linux systems. The current version was tested under
CentoOS 7.0.

One optional pre-requisite is DeconSeq (http://deconseq.sourceforge.net/),
which removes human reads. By default human reads are not removed.

## Installation

1. Create a working directory. Here we use ~/work as an example.
```
mkdir -p ~/work
cd ~/work
```

2. Get source code from github.
```
git clone https://github.com/eitanyaffe/hpipe.git
```

3. Set the environment variable HPIPE_DIR to point to the hpipe dir.
For example, if using bash, you can add the following line to your .bashrc:
```
export HPIPE_DIR=~/work/hpipe
```

4. Add the script `$HPIPE_DIR/hpipe` to your path. To do that either copy the
script to some common directory (e.g. /usr/local/bin), or add $HPIPE_DIR to your
path.

## Quick Start

Check hpipe on a test dataset:

1. Start an hpipe docker container.
```
hpipe start -c config/template/basic.cfg
```

2. Infer anchor-union pairs. This takes ~20 minutes on an 1TB/80-core machine.
```
hpipe run -c config/template/basic.cfg
```

3. When done close the container.
```
hpipe stop -c config/template/basic.cfg
```

Output files are generated under `$HPIPE_DIR/output/project1/result`.

## Input and Output

Input and output files and directories are defined using a single configuration
directory per project, see config/template for an example. To start a new project
you can copy the directory config/template and modify the various files there,
as explained below.

The configuration directory contains 3 files:
* project_id: Unique string identifier of project `PROJECT_ID`.
* path_vars: Paths to input and output directories.
* user_vars: Customizable parameters.
* basic.cfg: Experimental design. (Typically, do not edit this file.)

### Input

The pipeline uses as input two raw reads of DNA libraries.

* Assembly library, used to generate a metagenomic assembly. The location of the
paired fastq files is determined by the `ASSEMBLY_INPUTDIR` in the path_vars file.

* Hi-C library, used to infer the background model and infer anchor/genome
pairs. The location of the paired fastq files is determined by the `HIC_INPUTDIR` in
the path_vars file.

User-defined parameters in user_vars include:
* `CUTTER_TITLE`: The restriction enzyme used.
* `CUTTER_SITE`: The cutter site.
* `REMOVE_HUMAN`: [T|F] Should remove human reads using Deconseq.

### Output

The pipeline outputs all files under `BASE_OUTDIR`, defined in path_vars.
The output of each project is placed under `BASE_OUTDIR/PROJECT_ID/result`.
Temporary files are placed under `BASE_TMPDIR`.

Output files include:

* contig.fasta: fasta file of the metagenomic assembly.
* contig.table: table with all contigs and their length (in bp).
* anchor.table: table with all genome anchors, fields:
  * contig: contig identifier.
  * anchor: anchor identifier.
* contig_anchor.table: table with all genome unions, fields:
  * contig: contig identifier.
  * anchor: anchor identifier.
  * observed: number of contacts between contig and anchor.
  * expected: number of expected spurious contacts between contig and anchor.
  * score: enrichment of observed the expected, in log_10.
* model/fend.table: table with fragment ends, fields:
  * fend: fragment end identifier.
  * frag: fragment identifier.
  * strand: strand of fragment end.
  * contig: contig identifier.
  * coord: coordinate of fragment end.
  * frag_len: length of fragment.
  * abundance: abundance of contig.
  * anchor: anchor indentifier, if none equals 0.
  * abundance_bin: abundance bin.
  * frag_len_bin: fragment length bin.
* model/fragment_length.f, model/abundance.f: correction matrices.
* model/fragment_length.bins, model/abundance.bins: ranges of the bins.
* model/prior: constant prior probability value.
