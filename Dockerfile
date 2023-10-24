FROM continuumio/miniconda3

#PREAMBLE

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
    && apt-get --assume-yes upgrade

#MAIN

RUN conda env create -f cloud_environment.yaml \
    && conda clean --all --yes \
    && echo -e "#! /bin/bash\n\n# script to activate conda environment" > ~/.bashrc \
    && echo "export PS1='Docker> '" >> ~/.bashrc \
    && echo -e "\nconda activate TEsmall" >> ~/.bashrc

SHELL ["conda", "run", "-n", "TEsmall", "/bin/bash", "-c"]
RUN python setup.py install

#ENVIRONMENT

ENV LC_ALL C
ENV LANG C
ENV BASH_ENV ~/.bashrc