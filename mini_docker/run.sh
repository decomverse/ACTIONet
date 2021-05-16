#!/bin/bash
export PARENT_DOCKER=shmohammadi86/actionet:mini
export REGISTRY=shmohammadi86
export IMAGE_NAME=actionet
export VERSION=mini
docker run -it -e USER=ec2-user -e PASSWORD=insitro -e USERID=1000 -e GROUPID=1000 -e UMASK=022  -p 8888:8888  --rm -v /data:/data $PARENT_DOCKER


