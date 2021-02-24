PARENT_DOCKER=shmohammadi86/genomics:singlecell
REGISTRY=shmohammadi86
IMAGE_NAME=genomics
VERSION=actionet
all: build

build:
	docker build \
	  --tag $(IMAGE_NAME):${VERSION} --tag $(REGISTRY)/$(IMAGE_NAME):${VERSION} --build-arg PARENT_DOCKER=${PARENT_DOCKER} -t $(IMAGE_NAME) .
test:
	docker run -it -e USER=ec2-user -e PASSWORD=insitro -e USERID=1000 -e GROUPID=1000 -e UMASK=022 -p 8787:8787 -p 8888:8888  --rm  $(IMAGE_NAME):${VERSION}


