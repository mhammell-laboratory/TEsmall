FROM continuumio/miniconda3

#PREAMBLE

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
    && apt-get --assume-yes upgrade

#MAIN

RUN conda env create -f environment.yaml \
    && conda clean --all --yes

SHELL ["conda", "run", "-n", "TEsmall", "/bin/bash", "-c"]
RUN python setup.py install

RUN echo -e "source /opt/conda/etc/profile.d/conda.sh" > /home/genomics/setup.sh \
    && echo "conda activate TEsmall" >> /home/genomics/setup.sh

#ENVIRONMENT

ENV LC_ALL C
ENV LANG C
ENV BASH_ENV /home/genomics/setup.sh
