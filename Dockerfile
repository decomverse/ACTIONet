FROM ubuntu:18.04

RUN apt-get update && \
    apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get -yq install --no-install-recommends \
    wget ca-certificates libpq-dev gcc gcc-multilib g++ libsnappy-dev libhdf5-dev libsuitesparse-dev libnss3 xvfb make g++ libglpk40 git

COPY environment.yaml environment.yaml

# Install Miniconda with Python 3.9 into /opt
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# Enable Conda
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# Attach Conda to PATH
ENV PATH /opt/conda/bin:$PATH

# upgrade pip for faster dep resolution
RUN conda update -y -n base -c defaults conda

# Install ACTIONet environment
RUN conda env create -f environment.yaml && \
    echo "conda activate actionet" >> ~/.bashrc && \
    /opt/conda/bin/conda clean --all -fy && \
    rm -rf ~/.cache/pip

# Attach Conda to PATH
ENV PATH /opt/conda/envs/actionet/bin:$PATH

# Clean up after apt and conda
RUN apt-get clean && rm -rf /var/lib/apt/lists/*
RUN conda clean -tipy

# Set environment variables for Python
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8



RUN pip install jupyter && \
    jupyter notebook --generate-config -y

# Copy the entire repo & install ACTIONet
RUN mkdir ACTIONet
COPY . ACTIONet
WORKDIR ACTIONet
RUN git submodule update --init
RUN rm -rf build
RUN python setup.py build
RUN python setup.py develop

EXPOSE 8888
CMD screen -d -m bash -c "jupyter-lab --ip=0.0.0.0 --no-browser --allow-root --LabApp.token='' --LabApp.disable_check_xsrf=True"
