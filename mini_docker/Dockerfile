#FROM "shmohammadi86/genomics:singlecell"
FROM "jupyter/scipy-notebook"

USER root

#install python version
RUN apt-get update && apt-get install -yq --no-install-recommends libhdf5-dev libsuitesparse-dev libnss3 xvfb

#make sure cmake version is >=3.19
RUN conda install -c conda-forge cmake
RUN pip install git+https://github.com/shmohammadi86/ACTIONet@python-devel
EXPOSE 8888
CMD screen -d -m bash -c "jupyter-lab --ip=0.0.0.0 --no-browser --allow-root --LabApp.token='' --LabApp.disable_check_xsrf=True"
