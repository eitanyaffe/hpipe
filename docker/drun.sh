if [[ $# -eq 0 ]] ; then
    echo "usage: $(basename $0) <config_dir>"
    exit 0
fi
CFG_DIR=$1
# CFG_BASE=$(basename $CFG_DIR)

if [ ${CFG_DIR:0:1} != "/" ]; then
CFG_DIR=$PWD/$CFG_DIR
fi

LINK_FILE=$CFG_DIR/path_vars
PROJECT_ID=`cat ${CFG_DIR}/project_id`

# parse docker link file
echo "Reading volume file: $LINK_FILE"
echo "Providing directories in docker container:"
IO_PATHS="-v $HPIPE_DIR/pipeline:/hpipe -v $CFG_DIR:/hpipe/config"
while IFS='' read -r line || [[ -n "$line" ]]; do
    [ -z "$line" ] || [ "${line:0:1}" == "#" ] && continue
    array=($(eval echo ${line//=/ }))
    echo " /links/${array[0]} ==> ${array[1]}"
    IO_PATHS="$IO_PATHS -v $HOME:$HOME -v ${array[1]}:/links/${array[0]}"
    if [ ${array[0]} == "BASE_OUTDIR" ]; then 
	mkdir -p ${array[1]}
    fi
    if [ ${array[0]} == "BASE_TMPDIR" ]; then 
	mkdir -p ${array[1]}
    fi
    if [ ! -d ${array[1]} ]; then
	echo "Error: Directory ${array[1]} does not exist"
	exit 1
    fi
done <${LINK_FILE}

IMAGE_NAME=eitanyaffe/hpipe:v1.00
IMAGE_SHORT_NAME=hpipe
D_START_OPTS="-ti"
PERMISSIONS_PATHS="-v /etc/passwd:/etc/passwd -v /etc/shadow:/etc/shadow -v /etc/group:/etc/group"
CONTAINER_NAME=${IMAGE_SHORT_NAME}_${PROJECT_ID}_${USER}
HOST=${IMAGE_SHORT_NAME}_${PROJECT_ID}

CMD="docker run ${D_START_OPTS} -d -h ${HOST} -u ${USER} ${PERMISSIONS_PATHS} ${IO_PATHS} --name ${CONTAINER_NAME} ${IMAGE_NAME}"
echo "#" ${CMD}
${CMD}
