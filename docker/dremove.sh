if [[ $# -eq 0 ]] ; then
    echo "usage: $(basename $0) <config_dir>"
    exit 0
fi

CFG_DIR=$1
if [ ${CFG_DIR:0:1} != "/" ]; then
CFG_DIR=$PWD/$CFG_DIR
fi

PROJECT_ID=`cat ${CFG_DIR}/project_id`
IMAGE_NAME=hpipe
CONTAINER_NAME=${IMAGE_NAME}_${PROJECT_ID}_${USER}

CMD="docker stop ${CONTAINER_NAME}"
echo "#" ${CMD}
${CMD}

CMD="docker rm ${CONTAINER_NAME}"
echo "#" ${CMD}
${CMD}
