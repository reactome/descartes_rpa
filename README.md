[![Build Status](https://travis-ci.com/reactome/descartes_rpa.svg?branch=main)](https://travis-ci.com/reactome/descartes_rpa)
[![PyPI version](https://badge.fury.io/py/descartes-rpa.svg)](https://badge.fury.io/py/descartes-rpa)

# descartes_rpa

Python pipeline to extract pathway activity from single-cell clusters in a systematic manner, annotating each cluster with which pathway(s) it represents.

Pathways are annotated by extracting the differentialy expressed genes found in each clusters. Then, the pathways represented by these genes are retrieved using the Reactome analysis tools (explained [here](https://reactome.org/userguide/analysis)).

## Quickstart

[Here](https://nbviewer.jupyter.org/github/reactome/descartes_rpa/blob/main/demo/10x_data.ipynb) you can see examples of Single-Cell cluster pathway annotation using a 10x Genomics data set from the [main Scanpy tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html), or an [example](https://nbviewer.jupyter.org/github/reactome/descartes_rpa/blob/main/demo/mouse_data.ipynb) on Mouse Single-Cell data from the [Scanpy hematopoiesis in mouse tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html).

You can also check examples using the Descartes data set from human Single-Cell [tissues](https://github.com/reactome/descartes_rpa/tree/main/demo/tissues), such as [Pancreas](https://nbviewer.jupyter.org/github/reactome/descartes_rpa/blob/main/demo/tissues/Pancreas.ipynb) and [Liver](https://nbviewer.jupyter.org/github/reactome/descartes_rpa/blob/main/demo/tissues/Liver.ipynb).

## Installing

### Locally
Install the module using pip
```bash
pip install descartes_rpa
```

### Docker
Build image

```bash
docker-compose build descartes_rpa
```

Run the image

```bash
docker-compose run --rm descartes_rpa
```

## Descartes pathway data

Single-Cell cluster pathways data for each tissue in the descartes [Human Gene Expression During Development](https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/) database can be found in **[this link](https://drive.google.com/drive/u/2/folders/1TgcLzB5owOY8LsDkUqINLfXnTQ4a7t4N)**.
