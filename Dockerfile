FROM ubuntu:18.04
RUN apt-get update

RUN apt-get install -y wget pkg-config build-essential manpages-dev && rm -rf /var/lib/apt/lists/*
RUN gcc --version
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-aarch64.sh -b \
    && rm -f Miniconda3-latest-Linux-aarch64.sh 

COPY . /descartes_rpa/
WORKDIR /descartes_rpa/

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN conda --version
RUN conda env update -n base -f env.yml --prune

RUN python setup.py install

ENTRYPOINT ["/bin/bash"]
