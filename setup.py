from setuptools import setup


setup(
    name="descartes_rpa",
    version="1.0",
    description="descartes-rpa: Extract pathway features from Single-Cell",
    url="https://github.com/reactome/descartes-rpa",
    python_requires=">=3.9",
    author="Joao Luiz de Meirelles",
    author_email="jldemeirelles@gmail.com",
    packages=[
        "descartes_rpa",
        "descartes_rpa.convert",
    ],
    scripts=["descartes_rpa/convert/rds_to_mtx.R"],
    zip_safe=False,
    include_package_data=True
)
