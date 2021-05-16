FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive \
        LANG=C.UTF-8 LC_ALL=C.UTF-8 \
        PATH=/opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -yq wget curl screen unzip  bzip2 ca-certificates libglib2.0-0 git cmake build-essential software-properties-common libhdf5-dev libsuitesparse-dev libnss3 xvfb libcairo2-dev && \
    add-apt-repository ppa:cran/libgit2 && \
    apt-get update && \
    apt-get clean && \
    apt autoremove -y

RUN  ln -s /bin/tar /bin/gtar

#download the latest anaconda 
COPY files/base_environment.yml /tmp/base_environment.yml
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    conda config --system --set auto_update_conda false && \
    conda config --system --set show_channel_urls true && \
    conda update --quiet --yes conda && \
    conda env update -f /tmp/base_environment.yml -n base && \
    jupyter notebook --generate-config && \
    conda clean --all -ftipsy && \
    rm -rf ~/.cache/pip
    
##Copy config file
COPY files/jupyter_notebook_config.py /home/$DOCKER_USER/.jupyter/jupyter_notebook_config.py
COPY files/jupyter_notebook_config.py /root/.jupyter/jupyter_notebook_config.py

# We aren't running a GUI, so force matplotlib to use
# the non-interactive "Agg" backend for graphics.
# Run matplotlib once to build the font cache.
ENV MATPLOTLIBRC=/home/$DOCKER_USER/.config/matplotlib/matplotlibrc
RUN mkdir -p /home/$DOCKER_USER/.config/matplotlib && \
    echo "backend      : Agg" > /home/$DOCKER_USER/.config/matplotlib/matplotlibrc && \
    python -c "import matplotlib.pyplot"

# Instruct joblib to use disk for temporary files. Joblib defaults to
# /shm when that directory is present. In the Docker container, /shm is
# present but defaults to 64 MB.
# https://github.com/joblib/joblib/blob/0.11/joblib/parallel.py#L328L342
ENV JOBLIB_TEMP_FOLDER=/tmp
ENV PATH=/opt/conda/bin:$PATH \
    CONDARC=/opt/conda/.condarc \
    CONDA_DIR=/opt/conda

#install python actionet 
RUN pip install git+https://github.com/shmohammadi86/ACTIONet@python-release


#create an R conda environment
COPY files/r_environment.yml /tmp
RUN conda env update -f /tmp/r_environment.yml --name r_env && \    
    conda clean --all -f -y && \
    rm -rf ~/.cache/pip


ENV R_HOME=/opt/conda/envs/r_env/lib/R \
	CRAN=https://cloud.r-project.org

RUN echo "options(repos = c(CRAN = '${CRAN}'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site &&\
	echo  'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' >> ${R_HOME}/etc/Rprofile.site && \
	mkdir -p ${R_HOME}/site-library && \
	chown root:staff ${R_HOME}/site-library && \
	chmod g+ws ${R_HOME}/site-library && \
	echo "R_LIBS=\${R_LIBS-'${R_HOME}/site-library:${R_HOME}/library'}" >> ${R_HOME}/etc/Renviron && \
	echo "export MKL_NUM_THREADS=\`nproc\`" >> ~/.bashrc && \
	echo "MKL_THREADING_LAYER=GNU" >> ~/.bashrc


COPY files/*R /tmp/
#install R kernel 
RUN conda run -n r_env Rscript /tmp/install_irkernel.R
#install R actionet
RUN conda run -n r_env Rscript /tmp/install_actionet_deps.R
RUN conda run -n r_env Rscript /tmp/install_actionet.R 

EXPOSE 1-65535
CMD jupyter-lab --ip=0.0.0.0 --no-browser --allow-root --LabApp.token=''

