DOCKER_IMAGE='quay.io/pypa/manylinux_2_35_x86_64'
PLAT='manylinux_2_35_x86_64'
docker pull "$DOCKER_IMAGE"
docker container run -it --rm \
       -e PLAT=$PLAT \
       -v "$(pwd)":/io \
       "$DOCKER_IMAGE" /bin/bash 
