FROM continuumio/miniconda3

COPY . /descartes-rpa/
WORKDIR /descartes-rpa/

RUN conda env update -n base -f env.yml --prune

RUN python setup.py install

ENTRYPOINT ["/bin/bash"]
