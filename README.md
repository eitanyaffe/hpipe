# hpipe: A tool for the analysis of metagenomic Hi-C data

The hpipe pipeline processes metagenomic Hi-C maps. The core functionality is 
the inference of genome anchor/union pairs. Each genome anchor is a collection
contigs that are contained in the genomes of one ore more related strains, 
while each matching genome union is the combined genomes of the related strains. 
 
The pipeline was developed by Eitan Yaffe, at David Relman's lab, Stanford. 
It is distributed under the GNU General Public License v3.0. If you have
questions or comments please contact Eitan Yaffe at eitan.yaffe@gmail.com.

# Algorithm Overview

The hpipe pipeline steps include:

1. Library preprocess (Process raw reads).
  - Trim for quality
  - Remove adapters
  - Remove human reads (optional, requires Deconseq)

2. Assembly. 
  - Currently supports the minia assembler.
  - Generate contigs.

3. Anchor/Union infererence.
  - Computes seed anchors.
  - Infers background model over seed anchors.
  - Removes multi-anchor contigs until convergence.
  - Infers final model over anchors.
  - Extends each anchor into a matching genome union. 

# Prerequisites

The tool requires docker (https://www.docker.com) and a perl interpreter. All
work is done with a docker container, making the tool compatible with any 
system that supports docker.

Deconseq. 
TBD: Make it optional.

Genebank. For downloading ref genomes. 
TBD: Have local copy of ref genomes.

# Installation

Use a working directory of your choice. Here we use /work/hpipe as an example.

1. Get source code from github. 
%> mkdir -p /work
%> cd /work
%> git clone https://github.com/eitanyaffe/hpipe.git

This will download the hpile files into /work/hpipe. 

2. Set the environment variable HPIPE_DIR. For example, if using bash, add 
the following to your .bashrc:
export HPIPE_DIR=/work/hpipe

3. Add the $HPIPE_DIR/hpipe.pl wrapper script to your path. Either copy the
script to a common directory (e.g. /usr/local/bin) or add $HPIPE_DIR to your
path.

# Quick Start

To verify hpipe has been installed successfully you can run the following
commands.

1. start an hpipe container
%> ./hpipe.pl start -pdir pipeline -cdir config/ref -cfg config/ref/n5.cfg

2. generate simulated shotgun and Hi-C reads 
%> ./hpipe.pl run -pdir pipeline -cdir config/ref -cfg config/ref/n5.cfg -step \
   pp_simulate

3. infer anchor-union pairs
%> ./hpipe.pl run -pdir pipeline -cdir config/ref -cfg config/ref/n5.cfg -step \
   pp_basic

4. plot all figures
%> ./hpipe.pl run -pdir pipeline -cdir config/ref -cfg config/ref/n5.cfg -step \
   pl_basic

# Usage

The various pipeline parameters, including paths and values, are adjusted via
a configuration file. 

A) Edit the configuration file.

B) Starting/Stopping a container.

C) Running commands.

# FAQ

Below are solutions for common problems you may encounter.

--------------------------------------------------------------------------------

Q1: I am getting the following error message:
"Error response from daemon: Conflict. The name "hpipe_ref_eitany" is already 
in use by container fc86bf399649. You have to delete (or rename) that container 
to be able to reuse that name."

S1: A docker container was not left hanging for some reason. To remove the 
container run:
%> docker rm -f fc86bf399649

--------------------------------------------------------------------------------
