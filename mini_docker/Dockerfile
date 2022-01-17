FROM rocker/tidyverse
ARG GITHUB_PAT
LABEL name="actionet/tidyverse" \
	version=1.0.0 \
	vendor="ACTIONet" \
	maintainer="shahin.mohammadi@gmail.com" \
	description="ACTIONet single-cell framework." \
	license="GPLv2"


ENV DEBIAN_FRONTEND=noninteractive \
	APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 \
	R_VERSION=$R_VERSION

RUN apt-get update \
  && apt-get install -y --no-install-recommends libhdf5-dev libsuitesparse-dev libnss3 xvfb libglpk40
COPY files/* /tmp/
RUN /tmp/install_actionet.sh

CMD ["/init"]
