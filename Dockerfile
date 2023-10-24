FROM continuumio/miniconda3

#PREAMBLE

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
    && apt-get --assume-yes upgrade

#MAIN

RUN conda env create -f cloud_environment.yaml \
    && conda clean --all --yes

SHELL ["conda", "run", "-n", "TEsmall", "/bin/bash", "-c"]
RUN python setup.py install \
    && rm -rf /var/cache/apk/* \
    && rm -rf /tmp/*

#ENVIRONMENT

ENV LC_ALL C
ENV LANG C
ENV PATH /opt/conda/envs/TEsmall/bin:$PATH

