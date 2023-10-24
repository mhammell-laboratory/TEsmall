FROM continuumio/miniconda3

#PREAMBLE

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
    && apt-get --assume-yes upgrade

#MAIN

RUN conda env update -n base -f cloud_environment.yaml \
    && conda clean --all --yes \
    && python setup.py install

#ENVIRONMENT

ENV LC_ALL C
ENV LANG C
