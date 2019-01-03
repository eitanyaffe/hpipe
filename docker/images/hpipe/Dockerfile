FROM eitanyaffe/hpipe-base:v1.1

ENV MAKESHIFT_ROOT /ms_root

RUN mkdir ${MAKESHIFT_ROOT}
WORKDIR ${MAKESHIFT_ROOT}

# clone code from git
RUN \
git clone https://github.com/eitanyaffe/makeshift-core.git && \
git clone https://github.com/eitanyaffe/makeshift-modules.git && \
git clone https://github.com/eitanyaffe/hpipe

# compile binaries
WORKDIR ${MAKESHIFT_ROOT}/hpipe/pipeline
RUN make init

