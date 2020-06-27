ARG base_image=shmohammadi86/r4mkl_bioc:latest
FROM "$base_image"

LABEL name="shmohammadi86/actionet:latest" \
      version=1.0.0 \
      vendor="ACTIONet" \
      maintainer="shahin.mohammadi@gmail.com" \
      description="ACTIONet single-cell framework." \
      license="MIT"
      

ENV DEBIAN_FRONTEND noninteractive
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1

RUN Rscript -e "devtools::install_github('shmohammadi86/ACTIONet', ref = 'R-release')":

# Init command for s6-overlay
CMD ["/init"]
