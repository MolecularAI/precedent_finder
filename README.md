# Precedent finder

Example code for the Precedent finder project. The code works with the USPTO dataset.

The code was made for providing an open-source version of the tool, it was never intended to be production ready. Nor will it be maintained.


## Prerequisites

Before you begin, ensure you have met the following requirements:

* Linux, Windows or macOS platforms are supported - as long as the dependencies are supported on these platforms.

* You have installed [anaconda](https://www.anaconda.com/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) with python 3.8-3.11

The tool has been developed on a Linux platform.


## Installation

First clone the repository using Git.

Then execute the following commands in the root of the repository

    conda env create -f env.yml
    conda activate pfenv

Finally, you need to download the pre-processed USPTO dataset to the `data` folder. The files are located at [Zenodo](https://zenodo.org/records/15798247) and can be downloaded with

    python data/download_data.py


## Usage

To use the Precedent Finder tool use the `precedent_finder.py` script

    conda activate pfenv
    python precedent_finder.py --smiles REACTION_SMILES --output precedents.csv

or in interactive mode

    conda activate pfenv
    python precedent_finder.py

in this mode, you will be ask for one or more reaction SMILES and the results are saved to individual CSV-files.

**Note: the reaction SMILES need to be atom-mapped.** This can accomplished with the [rxnmapper](https://github.com/rxn4chemistry/rxnmapper/) project.

The `examples` folder contains worked examples with a Jupyter notebook showing how the output can be analyzed/visualized.


## Pre-processing data

To re-produce the pre-processing of the USPTO data, you can to do the following

First change to the `data` folder

    cd data

Second, perform the pipelines in the `rxnutils` package to process the USPTO data set, as described [here](https://molecularai.github.io/reaction_utils/uspto.html)

Finally, run the pipeline in the `data` folder with something like this

    conda activate pfenv
    python setup_pipeline.py run --nbatches 200  --max-workers 8 --max-num-splits 200


## Contributors

* [Samuel Genheden](https://www.github.com/SGenheden)
* [Christoph Bauer](https://www.github.com/CBA087)
* Thierry Kogej
* Per-Ola Norrby

The contributors have limited time for support questions.

## License

The software is licensed under the Apache 2.0 license (see LICENSE file), and is free and provided as-is.


## References