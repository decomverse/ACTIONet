FROM shmohammadi86/r-bioc

LABEL name="shmohammadi86/actionet:latest" \
      version=1.0.0 \
      vendor="ACTIONet" \
      maintainer="shahin.mohammadi@gmail.com" \
      description="ACTIONet single-cell framework." \
      license="GPLv2"
      
ENV DEBIAN_FRONTEND noninteractive

# Install ACTIONet
RUN Rscript -e 'update.packages(ask = FALSE) \
				devtools::install_github("shmohammadi86/SCINET", ref = "master") \
				devtools::install_github("shmohammadi86/NetLibR", ref = "master") \
				devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")'


CMD ["/init"]
