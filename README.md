# hpipe: A tool for the analysis of metagenomic Hi-C data

The hpipe tool processes metagenomic Hi-C maps to infer properties of the
underlying genomes, in the form of genome anchor/union pairs. Each genome 
anchor is a collection contigs that are contained in the genomes of one 
or more related strains, while each matching genome union is the combined 
genomes of the strains. 
 
The pipeline was developed by Eitan Yaffe, at David Relman's lab, Stanford. 
It is distributed under the GNU General Public License v3.0. If you have
questions or comments please contact Eitan Yaffe at eitan.yaffe@gmail.com.

--------------------------------------------------------------------------------
## Algorithm Overview

The major pipeline steps include:

1. Library preprocess.
  * Trim raw reads for quality.
  * Remove sequencing adapters.
  * Remove human reads (optional, requires a locally installed Deconseq pipeline).

2. Assembly. 
  * Currently supports the minia assembler.
  * Generate contigs.

3. Anchor/Union infererence.
  * Computes seed anchors.
  * Infers background model over seed anchors.
  * Removes multi-anchor contigs until convergence.
  * Infers final model over anchors.
  * Extends each anchor into a matching genome union. 

--------------------------------------------------------------------------------
## Prerequisites

The hpipe tool requires docker (https://www.docker.com) and a perl interpreter. 
All work is done within a docker container, making the tool compatible with most 
linux systems. The tool was tested under CentoOS 6.9.

Deconseq is a pipeline for the removal of human reads. It can be downloaded 
from http://deconseq.sourceforge.net/. The use of Deconseq is optional.

--------------------------------------------------------------------------------
## Installation

Select a working directory. Here we use /work/hpipe as an example.

1. Get source code from github. 
```
%> mkdir -p /work
%> cd /work
%> git clone https://github.com/eitanyaffe/hpipe.git
```

2. Set the environment variable HPIPE_DIR. For example, if using bash, add 
the following line to your .bashrc:
```
export HPIPE_DIR=/work/hpipe
```

3. Add the `$HPIPE_DIR/hpipe.pl` wrapper script to your path. Either copy the
script to a common directory (e.g. /usr/local/bin) or add $HPIPE_DIR to your
path.

--------------------------------------------------------------------------------
## Input and Output

### Input

The pipeline gets as input two DNA libraries. 
* Assembly library, used to generate a metagenomic assembly.
* Hi-C library, used to infer the background model and infer anchor/genome 
pairs.

--------------------------------------------------------------------------------
## Quick Start

To verify hpipe has been successfully installed run the following:

1. Start an hpipe container
```
%> ./hpipe.pl start -c config/example/example.cfg
```

2. Infer anchor-union pairs. This should take around 10 minutes.
```
%> ./hpipe.pl step -c config/example/example.cfg -s pp_basic
```

3. Close the container.
```
%> ./hpipe.pl stop -c config/example/example.cfg
```

All output files are generated under projects/project1.

--------------------------------------------------------------------------------
## Usage

The various pipeline parameters, including paths and values, are adjusted via
a configuration file. 

A) Prepare a configuration directory.

B) Managing containers.

C) Running commands.

--------------------------------------------------------------------------------
## Examples

### Beitel dataset

The Beitel dataset ([Beitel et al., PeerJ 2014](https://peerj.com/articles/415/)) contains
a Hi-C map for a synthetic mixture of 5 genomes. Since reference genomes are 
available the dataset comes without an assembly library. Therefore, to generate
an assembly, we simulate an assembly library from the reference genomes.

First download the Hi-C map from Sequence Read Archives (SRX377733).

To simulate the standard shotgun reads run:
```
%> ./hpipe.pl -c config/beitel/beitel.cfg -s pp_simulate
```

Then proceed normally to infer anchor/union pairs:
```
%> ./hpipe.pl -c config/beitel/beitel.cfg -s pp_basic
```

Finally you can plot the results:
```
%> ./hpipe.pl -c config/beitel/beitel.cfg -s pl_plot_basic
```

--------------------------------------------------------------------------------
## FAQ

Below are solutions to common problems.

--------------------------------------------------------------------------------


Q1: I am getting the following error message:
"Error response from daemon: Conflict. The name "hpipe_ref_eitany" is already 
in use by container fc86bf399649. You have to delete (or rename) that container 
to be able to reuse that name."

S1: A docker container was not left hanging for some reason. To remove the 
container run:
```
%> docker rm -f fc86bf399649
```

--------------------------------------------------------------------------------

