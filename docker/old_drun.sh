if [[ $# -eq 0 ]] ; then
    echo 'usage <config_file> <output dir> <output tmp dir> <target>'
    exit 0
fi
CONFIG_FILE=$1
OUTDIR=$2
TMPDIR=$3
TARGET=$4

CONFIG_DIR=$PWD/$(dirname $CONFIG_FILE)
PATHS_FILE=$CONFIG_DIR/docker_paths.cfg

echo "config file: $CONFIG_FILE"
echo "volume file: $PATHS_FILE"
echo "outdir: $OUTDIR"
echo "tmpdir: $TMPDIR"
echo "target: $TARGET"

IMAGE_NAME=hpipe_base
D_RUN_OPTS="-it"
PERMISSIONS_PATHS="-v /etc/passwd:/etc/passwd -v /etc/shadow:/etc/shadow -v /etc/group:/etc/group"
CONTAINER_NAME=${USER}_${IMAGE_NAME}

IO_PATHS="-v $CONFIG_DIR:/hpipe/config -v $OUTDIR:/hpipe/output -v $TMPDIR:/hpipe/tmp"
. $PATHS_FILE

MAKE_OPTS="-n"

echo "container: ${CONTAINER_NAME}"
if docker ps -a -f name=${CONTAINER_NAME} | grep -q ${CONTAINER_NAME}
then
    echo "removing container: ${CONTAINER_NAME}"
    # docker rm -f ${CONTAINER_NAME}
fi
# CMD="docker run ${D_RUN_OPTS} --rm -u ${USER} ${PERMISSIONS_PATHS} ${IO_PATHS} --name ${CONTAINER_NAME} ${IMAGE_NAME} make $TARGET $MAKE_OPTS"
CMD="docker run ${D_RUN_OPTS} --rm -u ${USER} ${PERMISSIONS_PATHS} ${IO_PATHS} --name ${CONTAINER_NAME} ${IMAGE_NAME} bash"
echo ${CMD}
${CMD}
