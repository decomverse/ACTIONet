FROM rocker/tidyverse

LABEL name="actionet/tidyverse" \
	version=1.0.0 \
	vendor="ACTIONet" \
	maintainer="shahin.mohammadi@gmail.com" \
	description="ACTIONet single-cell framework." \
	license="GPLv2"


ENV DEBIAN_FRONTEND=noninteractive \
	APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 \
	R_VERSION=$R_VERSION

COPY files/* /tmp/
RUN /tmp/install_actionet.sh

CMD ["/init"]
