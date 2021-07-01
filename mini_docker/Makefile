PARENT_DOCKER=actionet/actionet:mini
REGISTRY=actionet
IMAGE_NAME=actionet
VERSION=mini
all: build

build:
	docker build \
	  --tag $(IMAGE_NAME):${VERSION} --tag $(REGISTRY)/$(IMAGE_NAME):${VERSION} -t $(IMAGE_NAME) .
test:
	docker run -it -e USER=ec2-user -e PASSWORD=insitro -e USERID=1000 -e GROUPID=1000 -e UMASK=022 -p 8888:8888  --rm  $(IMAGE_NAME):${VERSION}

