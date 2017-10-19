if [[ $# -eq 0 ]] ; then
    echo "usage: $(basename $0) <hpipe_dir> <config_dir>"
    exit 0
fi
HPIPE_DIR=$1
CFG_DIR=$2

# add current dir if not absolute
if [ ${HPIPE_DIR:0:1} != "/" ]; then
HPIPE_DIR=$PWD/$HPIPE_DIR
fi
if [ ${CFG_DIR:0:1} != "/" ]; then
CFG_DIR=$PWD/$CFG_DIR
fi

LINK_FILE=$CFG_DIR/path_vars
PROJECT_ID=`cat ${CFG_DIR}/project_id`

# parse docker link file
echo "volume file: $LINK_FILE"
IO_PATHS="-v $HPIPE_DIR:/hpipe"
while read p; do
    [ -z "$p" ] || [ "${p:0:1}" == "#" ] && continue
    echo "line: $p"
    array=(${p//=/ })
    IO_PATHS="$IO_PATHS -v $HOME:$HOME -v ${array[1]}:/links/${array[0]}"
done <${LINK_FILE}

IMAGE_NAME=hpipe
D_START_OPTS="-ti --rm"
PERMISSIONS_PATHS="-v /etc/passwd:/etc/passwd -v /etc/shadow:/etc/shadow -v /etc/group:/etc/group"
CONTAINER_NAME=${IMAGE_NAME}_${PROJECT_ID}_${USER}
HOST=${IMAGE_NAME}_${PROJECT_ID}

CMD="docker run ${D_START_OPTS} -h ${HOST} -u ${USER} ${PERMISSIONS_PATHS} ${IO_PATHS} --name ${CONTAINER_NAME} ${IMAGE_NAME}"
echo "#" ${CMD}
${CMD}

