FROM shmohammadi86/r-bioc

LABEL name="shmohammadi86/actionet:latest" \
      version=1.0.0 \
      vendor="ACTIONet" \
      maintainer="shahin.mohammadi@gmail.com" \
      description="ACTIONet single-cell framework." \
      license="GPLv2"
      
ENV DEBIAN_FRONTEND noninteractive

# Install ACTIONet
RUN Rscript -e 'devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")'

## Install nextflow
RUN apt-get update && \
    apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get -y install --no-install-recommends \
        openjdk-8-jre-headless openssh-client cmake && \
    rm -rf /var/lib/apt/lists/*
RUN cd /tmp && export NXF_VER=19.07.0 && curl https://get.nextflow.io | bash && mv nextflow /usr/local/bin

## Install AWS CLI

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
  unzip awscliv2.zip && \
  sudo ./aws/install

CMD ["/init"]
