FROM centos:7

# R
RUN \
  yum update -y && \
  yum upgrade -y && \
  yum install -y epel-release && \
  yum install -y R

# various tools
RUN \
yum install -y http://mirror.centos.org/centos/7/os/x86_64/Packages/time-1.7-45.el7.x86_64.rpm && \
yum install -y wget && \
yum install -y git && \
yum install -y parallel && echo 'will cite' | parallel --bibtex && \
yum -y install emacs-nox && \
yum install -y perl-CPAN && \
yum install -y perl-Compress-Raw-Zlib && \
curl -L http://cpanmin.us | perl - App::cpanminus && \
cpanm Digest::MD5 && cpanm Graph && cpanm Switch

########################################################################
# SeqPrep
########################################################################

RUN git clone https://github.com/jstjohn/SeqPrep.git /tools/SeqPrep && make -C /tools/SeqPrep
ENV SEQPREP /tools/SeqPrep/SeqPrep

########################################################################
# sickle
########################################################################

RUN git clone https://github.com/najoshi/sickle.git /tools/sickle && make -C /tools/sickle
ENV SICKLE /tools/sickle/sickle

########################################################################
# MUMmer
########################################################################

ENV MUMMER_VERSION 3.23
ENV MUMMER_DIR /tools/mummer/MUMmer${MUMMER_VERSION}
WORKDIR /tools/mummer
RUN wget https://downloads.sourceforge.net/project/mummer/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz &&  \
tar xvf MUMmer${MUMMER_VERSION}.tar.gz
WORKDIR /tools/mummer/MUMmer${MUMMER_VERSION}
RUN make

########################################################################
# bwa
########################################################################

ENV BWA_VERSION 0.7.12-1.el7
RUN yum install -y bwa-${BWA_VERSION}
ENV BWA_BIN bwa

########################################################################
# diamond
########################################################################

ENV DIAMOND_VERSION v0.9.10
WORKDIR /tools/diamond
RUN wget http://github.com/bbuchfink/diamond/releases/download/${DIAMOND_VERSION}/diamond-linux64.tar.gz && \
tar xvf diamond-linux64.tar.gz
ENV DIAMOND_BIN /tools/diamond/diamond

########################################################################
# megahit
########################################################################

ENV MEGAHIT_VERSION v1.1.4
ENV MEGAHIT_DIR /tools/megahit
WORKDIR ${MEGAHIT_DIR}
RUN git clone -b ${MEGAHIT_VERSION} https://github.com/voutcn/megahit.git ${MEGAHIT_DIR}
RUN yum install -y libstdc++-static
RUN make
ENV MEGAHIT ${MEGAHIT_DIR}/megahit

########################################################################
# R packages
########################################################################

WORKDIR /tools
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > /tools/.Rprofile && \
Rscript -e "install.packages(c('gplots', 'R.utils', 'igraph', 'maxLik', 'fastcluster', 'mixtools'))"

