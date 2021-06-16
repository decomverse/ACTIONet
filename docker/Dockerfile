ARG PARENT_DOCKER
FROM "$PARENT_DOCKER"
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS="true"
ENV GITHUB_PAT=

#install R version 
COPY files/* /tmp/
RUN Rscript /tmp/install_actionet.R

#install python version
RUN apt-get install -y libhdf5-dev libsuitesparse-dev libnss3 xvfb

#make sure cmake version is >=3.19
RUN conda install -c conda-forge cmake
RUN pip install git+https://github.com/shmohammadi86/ACTIONet@python-release
EXPOSE 8787 8888
CMD screen -d -m bash -c "jupyter-lab --ip=0.0.0.0 --no-browser --allow-root --LabApp.token=''" && /init
