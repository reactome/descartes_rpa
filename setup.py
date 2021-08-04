from setuptools import setup


with open("README.md", "r") as f:
    long_description = f.read()


setup(
    name="descartes_rpa",
    version="1.2",
    description="descartes_rpa: Extract pathway features from Single-Cell",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/reactome/descartes_rpa",
    python_requires=">=3.8",
    author="Joao Luiz de Meirelles",
    author_email="jldemeirelles@gmail.com",
    packages=[
        "descartes_rpa",
        "descartes_rpa.convert",
        "descartes_rpa.fetch",
        "descartes_rpa.io",
        "descartes_rpa.pl",
        "descartes_rpa.analyze",
        "descartes_rpa.test"
    ],
    install_requires=[
        "scanpy==1.7.2",
        "requests==2.25.1",
        "aiohttp==3.7.4",
        "flake8==3.9.2",
        "loompy==3.0.6",
        "reactome2py==3.0.0",
        "seaborn==0.11.1",
        "scikit-learn==0.24.2",
        "statsmodels==0.12.2",
        "numba==0.53.1",
        "tables==3.6.1",
        "python-igraph==0.9.1",
        "leidenalg==0.8.4",
        "UpSetPlot==0.5.0"
    ],
    include_package_data=True
)
