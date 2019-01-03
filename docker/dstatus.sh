if [[ $# -eq 0 ]] ; then
    echo "usage: $(basename $0) <config_dir>"
    exit 0
fi
CFG_DIR=$1

if [ ! -d ${CFG_DIR} ]; then
    echo "Error: directory ${CFG_DIR} does not exist"
    exit 1
fi

if [ ${CFG_DIR:0:1} != "/" ]; then
CFG_DIR=$PWD/$CFG_DIR
fi

PROJECT_ID=`cat ${CFG_DIR}/project_id`
IMAGE_SHORT_NAME=hpipe
CONTAINER_NAME=${IMAGE_SHORT_NAME}_${PROJECT_ID}_${USER}

CMD="docker ps -f name=${CONTAINER_NAME} -q"
# echo "#" ${CMD}
ID="$($CMD)"
if [[ ${ID} != "" ]]; then
    echo "${CONTAINER_NAME} container is running, ID: ${ID}"
else
    echo "${CONTAINER_NAME} container is not running"
fi
