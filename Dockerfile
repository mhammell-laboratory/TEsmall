FROM continuumio/miniconda3

#PREAMBLE

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
               && apt-gel --assume-yes upgrade

#MAIN

RUN conda env create -f cloud_environment.yaml
RUN conda clean --all --yes
SHELL ["conda", "run", "-n", "TEsmall", "/bin/bash", "-c"]
RUN python setup.py install

