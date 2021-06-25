from setuptools import setup


with open("README.md", "r") as f:
    long_description = f.read()


setup(
    name="descartes_rpa",
    version="1.0",
    description="descartes-rpa: Extract pathway features from Single-Cell",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/reactome/descartes-rpa",
    python_requires=">=3.9",
    author="Joao Luiz de Meirelles",
    author_email="jldemeirelles@gmail.com",
    packages=[
        "descartes_rpa",
        "descartes_rpa.convert",
        "descartes_rpa.fetch",
        "descartes_rpa.analyze",
        "descartes_rpa.test"
    ],
    zip_safe=False,
    include_package_data=True
)
