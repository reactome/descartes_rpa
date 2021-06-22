FROM continuumio/miniconda3

COPY . /descartes_rpa/
WORKDIR /descartes_rpa/

RUN conda env update -n base -f env.yml --prune

RUN python setup.py install

ENTRYPOINT ["/bin/bash"]
