if [[ $# -eq 0 ]] ; then
    echo "usage: $(basename $0) <config_dir> <command>"
    exit 0
fi

CFG_DIR=$1
COMMAND=$2

# add current dir if not absolute
if [ ${CFG_DIR:0:1} != "/" ]; then
CFG_DIR=$PWD/$CFG_DIR
fi

PROJECT_ID=`cat ${CFG_DIR}/project_id`
IMAGE_NAME=hpipe
D_EXEC_OPTS="-ti"
CONTAINER_NAME=${IMAGE_NAME}_${PROJECT_ID}_${USER}

CMD="docker exec ${D_EXEC_OPTS} -u ${USER} ${CONTAINER_NAME} make $COMMAND"
echo "#" ${CMD}
${CMD}
