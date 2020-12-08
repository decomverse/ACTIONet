ARG base_image=nuest/mro:4.0.2-rstudio
FROM "$base_image"

LABEL name="shmohammadi86/actionet:latest" \
      version=1.0.0 \
      vendor="ACTIONet" \
      maintainer="shahin.mohammadi@gmail.com" \
      description="ACTIONet single-cell framework." \
      license="GPLv2"

RUN Rscript -e 'devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")'

CMD ["/init"]
