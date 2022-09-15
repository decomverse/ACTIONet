FROM "jupyter/scipy-notebook"

USER root

RUN apt-get update && apt-get install -yq --no-install-recommends libhdf5-dev libsuitesparse-dev libnss3 xvfb make g++ libglpk40

COPY environment.yaml environment.yaml

RUN conda config --system --set auto_update_conda false && \
    conda config --system --set show_channel_urls true && \
    conda update --quiet --yes conda && \
    conda env create -f environment.yml  && \
    jupyter notebook --generate-config && \
    conda clean --all -fy && \
    echo "conda activate actionet" >> ~/.bashrc && \
    rm -rf ~/.cache/pip
    
RUN python3 -m pip install types-setuptools
RUN mypy --install-types

# need to upgrade pip for faster dependency resolution
RUN pip install --upgrade pip

# Copy the entire repo & install ACTIONet
RUN mkdir ACTIONet
COPY . ACTIONet
RUN pip install -e ACTIONet

EXPOSE 8888
CMD screen -d -m bash -c "jupyter-lab --ip=0.0.0.0 --no-browser --allow-root --LabApp.token='' --LabApp.disable_check_xsrf=True"
