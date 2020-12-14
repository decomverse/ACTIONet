FROM nuest/mro:4.0.2-rstudio

LABEL name="shmohammadi86/actionet" \
      version=1.0.0 \
      vendor="ACTIONet" \
      maintainer="shahin.mohammadi@gmail.com" \
      description="ACTIONet single-cell framework." \
      license="GPLv2"

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	libhdf5-dev \
	libsuitesparse-dev \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

     
RUN mkdir ~/.R/ && echo "MAKEFLAGS+=-j \`nproc\`" >> ~/.R/Makevars && \
	echo "CURL_CA_BUNDLE=/opt/microsoft/ropen/$MRO_VERSION//lib64/R/lib/microsoft-r-cacert.pem" >> /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Renviron

RUN cp /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile.site /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile_old.site && head -n -7 /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile.site > /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile_new.site && sed -i 's/quiet <- /quiet <- TRUE #/' /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile_new.site && cp /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile_new.site /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile.site && echo "options("download.file.method" = \"libcurl\")" >> /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile.site

RUN R -e 'install.packages(c("devtools", "httr", "curl"))'
RUN R -e 'devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")'

RUN cp /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile_old.site /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile.site && rm /opt/microsoft/ropen/$MRO_VERSION/lib64/R/etc/Rprofile_old.site

CMD ["/init"]
